#!/usr/bin/env nextflow
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow run main.nf --tumor samples/ --normal samples/ --gc hg19.gc.wig --cen hg19.centromere.txt -profile docker
""".stripIndent()
}
params.genome = 'GRCh37'

// Read in params.genome files
ref_fasta = file(params.genomes[params.genome].ref_fasta, checkIfExists: true)
ref_index = Channel.fromPath(params.genomes[params.genome].ref_index, checkIfExists: true).collect()
gatk_mills = file(params.genomes[params.genome].gatk_mills, checkIfExists: true)
gatk_mills_index = file(params.genomes[params.genome].gatk_mills_index, checkIfExists: true)
gatk_1kg = file(params.genomes[params.genome].gatk_1kg, checkIfExists: true)
gatk_1kg_index = file(params.genomes[params.genome].gatk_1kg_index, checkIfExists: true)

// Read in command line input files
gc_window = file(params.gc_window, checkIfExists: true)
centromere_file = file(params.centromere, checkIfExists: true)
Channel.fromFilePairs(params.normal + '/*{1,2}.fastq.gz', flat: true, checkIfExists: true)
                       .map { tuple( it[0], "normal", it[1], it[2] ) }
                       .set {normal_samples}

Channel.fromFilePairs(params.tumor + '/*{1,2}.fastq.gz', flat: true, checkIfExists: true)
                       .map { tuple( it[0], "tumor", it[1], it[2] ) }
                       .set {tumor_samples}

fastqs = normal_samples.mix(tumor_samples)


process alignment {
    // Align fastqs, sort and index
    label 'alignment'

    tag "${sample_id}-${sample_type}"

    input:
        path ref_fasta
        path ref_index
        tuple val(sample_id), val(sample_type), file(fastq1), file(fastq2) from fastqs

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.raw.bam") into raw_bams

    // publishDir params.output, mode: 'copy', overwrite: true

    // -K process INT input bases in each batch regardless of nThreads (for reproducibility)
    script:
    """ 
    bwa mem \
       -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:NA\\tSM:${sample_id}\\t" \
       -K 100000000 \
       -t ${task.cpus}  \
       ${ref_fasta} ${fastq1} ${fastq2} 2> log.txt \
     | samtools sort -t${task.cpus} -m4G - -o ${sample_id}.${sample_type}.raw.bam
     """
}


process picard_remove_duplicates {
    label 'picard'

    tag "${sample_id}-${sample_type}"

    input:
        tuple val(sample_id), val(sample_type), file(bam_file) from raw_bams

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.rmdup.bam"), file("${sample_id}.${sample_type}.rmdup.bai") into rmdup_bams

    //publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ -Dpicard.useLegacyParser=false \
    MarkDuplicates \
    -INPUT ${bam_file} \
    -OUTPUT ${sample_id}.${sample_type}.rmdup.bam \
    -METRICS_FILE ${sample_id}.${sample_type}.quality_metrics \
    -REMOVE_DUPLICATES true \
    -ASSUME_SORTED true \
    -VALIDATION_STRINGENCY SILENT \
    -CREATE_INDEX true 2> picard_rmdupes.log
    """
}

process gatk_bqsr {
    label 'gatk'

    tag "${sample_id}-${sample_type}"

    input:
        path ref_fasta
        path ref_index
        path gatk_mills
        path gatk_mills_index
        path gatk_1kg
        path gatk_1kg_index
        tuple val(sample_id), val(sample_type), file(bam_file), file(bam_bai) from rmdup_bams

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.bqsr.bam") into bqsr_bams

    //publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    BaseRecalibrator \
    --reference ${ref_fasta} \
    --input ${bam_file} \
    --known-sites ${gatk_mills} \
    --known-sites ${gatk_1kg} \
    --output ${sample_id}.${sample_type}.recal_table

    gatk ApplyBQSR \
    --reference ${ref_fasta} \
    --input ${bam_file} \
    --bqsr-recal-file ${sample_id}.${sample_type}.recal_table \
    --output ${sample_id}.${sample_type}.bqsr.bam
    """
}


process samtools_final_bam {
    label 'samtools'

    tag "${sample_id}-${sample_type}"

    input:
        tuple val(sample_id), val(sample_type), file(bqsr_bam) from bqsr_bams

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.final.bam") into final_bams
        file("*.bai")

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    samtools sort ${bqsr_bam} -o ${sample_id}.${sample_type}.final.bam

    samtools index ${sample_id}.${sample_type}.final.bam
    """
}

process samtools_mpileup {
    label 'samtools'

    tag "${sample_id}-${sample_type}"


    input:
        path ref_fasta
        tuple val(sample_id), val(sample_type), file(final_bam) from final_bams

    output:
        tuple val(sample_id), file("${sample_id}.${sample_type}.mpileup") into mpileups

    publishDir params.output, mode: 'copy', overwrite: true

    // -l ${bed_file}
    // -B --no-BAQ
    // -E --redo-BAQ
    script:
    """
    samtools mpileup -f ${ref_fasta} -d 1000000 -A -B ${final_bam} > ${sample_id}.${sample_type}.mpileup
    """
}

process sequenza_pileup2seqz {
    label 'sequenza'

    tag "${sample_id}"

    echo true

    paired_mpileups = mpileups
        .groupTuple(size: 2)
        .map{
            sample_id, mpileup_files -> tuple( sample_id, mpileup_files.sort{ it.getName() } )
            }

    input:
        path gc_window
        tuple val(sample_id), file(mpileup_files) from paired_mpileups
        // [sample_id, [normal.mpileup, tumor.mpileup] ]

    output:
        tuple val(sample_id), file("${sample_id}.seqz") into sequenza_seqz

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    sequenza-utils bam2seqz -gc ${gc_window} -p -n ${mpileup_files[0]} -t ${mpileup_files[1]} -o ${sample_id}.seqz
    """
}

process sequenza_seqz_binning {
    label 'sequenza'

    tag "${sample_id}"

    //echo true

    input:
        tuple val(sample_id), file(seqz_gz) from sequenza_seqz.filter{ it[1].countLines() > 1 }
        // .filter{} -- seqz files that contain more than just a header line

    output:
        tuple val(sample_id), file("${sample_id}.binned.seqz.gz") into binned_seqz

    //publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    sequenza-utils seqz_binning -w 50 -s ${seqz_gz} -o ${sample_id}.binned.seqz.gz
    """
}

process sequenza_R {
    label 'sequenza'

    tag "${sample_id}"

    //echo true

    input:
        tuple val(sample_id), file(binned_seqz_gz) from binned_seqz

    output:
        tuple val(sample_id), file("${sample_id}.nitz.cellularity.txt"), file("${sample_id}.nitz.ploidy.txt"), file("${sample_id}.nitz.ave_depth.txt"), file("${sample_id}.nitz.copynumber_calls.txt"), file("${sample_id}_genome_view.pdf") into sequenza_R_files

    publishDir params.output, mode: 'copy', overwrite: true

    shell:
    '''
	#-----------------------------------------------------------------------------------------------
	#RUN SEQUENZA, R
	#-----------------------------------------------------------------------------------------------
	#Create an executable R script, run it and quit it!
	echo 'library(\"sequenza\")'>!{sample_id}.sequenza.r
	echo 'data.file <- \"!{binned_seqz_gz}\"' >> !{sample_id}.sequenza.r
	echo 'seqz.data <- read.seqz(data.file)' >> !{sample_id}.sequenza.r
	echo 'gc.stats <- gc.sample.stats(data.file)' >> !{sample_id}.sequenza.r
	echo 'test <- sequenza.extract(data.file)' >> !{sample_id}.sequenza.r
	echo 'CP.example <- sequenza.fit(test)' >> !{sample_id}.sequenza.r
	echo 'sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = \"!{sample_id}\", out.dir=\"./\")' >> !{sample_id}.sequenza.r
	echo 'cint <- get.ci(CP.example)' >> !{sample_id}.sequenza.r

	#Plot cellularity
	echo 'jpeg(\"!{sample_id}.nitz.cellularity.jpg\")' >> !{sample_id}.sequenza.r
	echo 'cp.plot(CP.example)' >> !{sample_id}.sequenza.r
	echo 'cp.plot.contours(CP.example, add = TRUE, likThresh=c(0.95))' >> !{sample_id}.sequenza.r
	echo 'dev.off()' >> !{sample_id}.sequenza.r

	#Call CNVs
	echo 'cellularity <- cint\$max.cellularity' >> !{sample_id}.sequenza.r
	echo 'ploidy <- cint\$max.ploidy' >> !{sample_id}.sequenza.r
    echo 'seg_table <- read.table(\"!{sample_id}_segments.txt\", header = TRUE, sep = \"\\t\", dec = \".\")' >> !{sample_id}.sequenza.r
	echo 'avg.depth.ratio <- mean(seg_table$depth.ratio)' >> !{sample_id}.sequenza.r

	#Save parameters to file
	echo 'cellularity' >> !{sample_id}.sequenza.r
	echo 'write(cellularity, file = \"!{sample_id}.nitz.cellularity.txt\")' >> !{sample_id}.sequenza.r
	echo 'write(ploidy, file = \"!{sample_id}.nitz.ploidy.txt\")' >>!{sample_id}.sequenza.r
	echo 'write(avg.depth.ratio, file = \"!{sample_id}.nitz.ave_depth.txt\")' >> !{sample_id}.sequenza.r

	#Detect variant alleles
	echo 'mut.tab <- na.exclude(do.call(rbind, test\$mutations))' >> !{sample_id}.sequenza.r
	echo 'mut.alleles <- mufreq.bayes(mufreq = mut.tab\$F, depth.ratio = mut.tab\$adjusted.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)' >> !{sample_id}.sequenza.r

	#Detect CN variation
	echo 'seg.tab <- na.exclude(do.call(rbind, test\$segments))' >> !{sample_id}.sequenza.r
	echo 'cn.alleles <- baf.bayes(Bf = seg.tab\$Bf, depth.ratio = seg.tab\$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)' >> !{sample_id}.sequenza.r
	echo 'seg.tab <- cbind(seg.tab, cn.alleles)' >>!{sample_id}.sequenza.r
	echo 'seg.tab' >> !{sample_id}.sequenza.r

	#write sequenza matrix to file, this will serve as input to loss score script's 2nd arg
	echo 'write.table(seg.tab, file = \"!{sample_id}.nitz.copynumber_calls.txt\", append = FALSE)' >> !{sample_id}.sequenza.r
	#exit
	echo 'q()' >> !{sample_id}.sequenza.r
	echo 'n' >> !{sample_id}.sequenza.r

	#execute the R script
	R --vanilla < !{sample_id}.sequenza.r
    '''
}

process loh_score {
    label 'sequenza'

    tag "${sample_id}"

    //echo true

    input:
        path centromere_file
        tuple val(sample_id), file(cellularity), file(ploidy), file(ave_depth), file(copynumber_calls), file(genome_pdf) from sequenza_R_files

    output:
        tuple val(sample_id), file("${sample_id}.nitz.score.txt") into scoring_output

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    LOH_score_chr_arms_V4.py ${centromere_file} ${copynumber_calls} ${sample_id}.nitz.score.txt 0.75
    echo "" >> ${sample_id}.nitz.score.txt
    echo -n "Estimated tumor cellularity: " >> ${sample_id}.nitz.score.txt
    cat ${cellularity} >> ${sample_id}.nitz.score.txt
    echo -n "Estimated ploidy: " >> ${sample_id}.nitz.score.txt
    cat ${ploidy} >> ${sample_id}.nitz.score.txt
    """
}
