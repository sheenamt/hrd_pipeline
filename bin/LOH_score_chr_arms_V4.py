#!/usr/bin/env python2
# coding: utf-8
# This version is an extention of Version 2 with 2 changes 1) analysis of X chromosome was removed 2) if cnt for a segment is  0, then the bases are not counted as lost bases
# This script does not consider cases where loss A or loss B over all available segments in a chromosomal arm/chromosome length span over 90% of the arm/chromosome. Such consistent behaviour over large areas of chromosome do NOT indicate HRD!
# If we refer to the Lancelet paper "Rucaparib in relapsed, platinum-sensitive high-grade ovarian carcinoma (ARIEL2 Part 1): an international, multicentre, open-label, phase 2 trial", its supplementary material states the following:
# To compute the percent genomic LOH for each tumour, LOH segments were inferred across the 22 autosomal chromosomes using the genome-wide aneuploidy/copy number profile and minor allele frequencies of the more than 3500 polymorphic SNPs sequenced in the Foundation Medicine’s NGS-based T5a assay.
# This profile was segmented and interpreted using allele frequencies of sequenced SNPs to estimate copy number (Ci) and minor allele count (Mi) at each segment (i). A segment was determined to have LOH if Ci ≠ 0 and Mi = 0.
# Two types of LOH segments were excluded from the calculation of percent genomic LOH: (1) LOH segments spanning ≥90% of a whole chromosome or chromosome arm, as these LOH events usually arise through non-HRD mechanisms (eg, mitotic nondisjunction6 ), and (2) regions in which LOH inference was ambiguous 

# With respect to the above assumptions, we do the following:
# 1. we first ensure that no segment in the copy number file spans over a centromere, if so, amend line
# 2. Only consider autosomal chormosome (No X chr)
# 3. Looking at each segment at a time, we 1) determine segment length 2) add segment length to chr_x_length (chrx can be chr1/chr2 etc) 3) determine if segment is in p/q arm
# 4. Depending on which chromosomal arm the segment is in, we 1) add segment length to chr_x_arm_length 
# 4. If the valid segment's count (A+B) is more than 0 (meaning this isn't a deletion/no coverage event) and either A or B is 0 (but not both), count the segment length toward 1) chr_x_arm_lossA/B 2) chr_x_lossA/B
# 5. Calculate the LOH score as follows:
#	For each chromosome: 
#	we add to total number of valid bases and total number of bases lost depending on the following conditions:
#	a) calculate A loss and B loss over entire chromosome and see if (lossA_chr_x AND lossB_chr_x) < 0.9  (equation used:  Aloss = chr_x_lossA/B / chr_x_length)
#		b) if calculate A loss and B loss over chromosomal arm and see if (lossA_chr_x_arm AND lossB_chr_x_arm) < 0.9 (equation used: loss_A_chr_x_arm = chr_x_arm_lossA/B / chr_x_arm_length)
#			c) if so, then add 1) arm length to total number of valid bases - denominator 2) total Aloss_arm and Bloss_arm to total number of bases lost - denominator 
#	and find the LOH score as total bases lost/ validBases (this is over ALL chromosomes)

# Note the chr_x_arm_length, chr_x_length are all total number of interrogatable bases through Sequenza's findings, not actual lengths of chromosomes/chromosomal arms.          

#This script does NOT consider chromosomes which have 90% or more loss over the entire chromosome or chromosomal arms that have 90% or more loss as these events are often not HRD

import re
import sys

centromere_info = {}
chromosome_info = {}

def centromereInfo(input_file):
	#Let us store the centromere start and stop positions
	with open(input_file, "r+") as centromere_file:
		header = centromere_file.readline()
		for line in centromere_file:
			chrom = line.split("\t")[1].replace("chr","")
			centromere_start = (line.split("\t")[2])
			centromere_end = (line.split("\t")[3])
			centromere_info[chrom] = [centromere_start, centromere_end]
		#print centromere_info



def checkCopyNoFile (input_file, modified_input_file):
	#the first thing to do is to go back to copynumber file and see if any segments cross over a centromere, if so, split the line
	data = []
	with open(input_file, "a+") as copyno:
		data.append(copyno.readline())
		for line in copyno:
			chromosome = (line.split(" ")[1].replace("\"","").replace("chr",""))
			start = line.split(" ")[2]
			stop = line.split(" ")[3]
	
			for chrom,positions in centromere_info.items():
				if chromosome == chrom:
					print "on chrom " + chromosome

					#Case1: if centromere is completely within segment
					if int(start) < int(positions[0]) and int(stop) > int(positions[1]):
						print "case1"
						data.append(line.replace(stop, positions[0]))
						data.append(line.replace(start, positions[1]))

					#Case2: if segment starts before the centromere but ends within the centromere
					elif  int(start) < int(positions[0]) and (int(stop) > int(positions[0]) and int(stop) <= int(positions[1])):
						print "case2"
						data.append(line.replace(stop, positions[0]))

					#Case3: if segment starts in the centromere but ends after the centromere
					elif  (int(start) >= int(positions[0]) and int(start) < int(positions[1])) and int(stop) > int(positions[1]):
						print "case3"
						data.append(line.replace(start, positions[1]))

					#Case4: if segment is completely within centromere
					elif int(start) >= int(positions[0]) and int(stop) <= int(positions[1]):
						print "case4"
						#do not append line

					#If none of these cases apply, just append line
					else:
						data.append(line)

	with open(modified_input_file, "w") as modified_file:
		modified_file.writelines(data)

						
validChromosomeArms = ","
#only consider autosomal chromosomes
chromosomesToTest = range(1,23)

total_bases = 0
total_lost_bases = 0

def calculateLOHscore(modified_file_path_file, chromosome, loss_threshold, p_arm_len = 0, p_arm_lossA = 0, p_arm_lossB = 0, q_arm_len = 0, q_arm_lossA = 0, q_arm_lossB = 0, chr_len = 0, chr_lossA = 0, chr_lossB = 0):
	print "====================================================================="
	print "processing chr: " + str(chromosome)
	global validChromosomeArms
	global total_bases
	global total_lost_bases
	
	with open(input_file, "r+") as seq_matrix:
		header = seq_matrix.readline()
		for line in seq_matrix:
			
			segmentLength = 0

			line =  re.sub("\"", "", line)
			line = re.sub("chr", "", line)
			line = line.replace("\n", "")

			outputArray = line.split(" ") #forms a list of strs

			if str(chromosome) == outputArray[1]:
				segmentLength = int(outputArray[3]) - int(outputArray[2])
				chr_len += segmentLength
				if int(outputArray[2]) < int(centromere_info[str(chromosome)][0]) and int(outputArray[3]) <= int(centromere_info[str(chromosome)][0]): #if segment start < chrom's centromere start and segment stop <= same
					p_arm_len += segmentLength #means segment in in upper arm of chromosome

					if int(outputArray[10]) > 0 : #cnt is more than 0 (if cnt is 0, this is either a deletion or no coverage, either way, it is not indicative of LOH event)
						if outputArray[11] == "0": #estimated number of A alleles is 0
							p_arm_lossA += segmentLength
							chr_lossA += segmentLength
						if outputArray[12] == "0": #estimated numbr of B alleles is 0
							p_arm_lossB += segmentLength
							chr_lossB += segmentLength

				elif int(outputArray[2]) >= int(centromere_info[str(chromosome)][1]) and int(outputArray[3]) > int(centromere_info[str(chromosome)][1]): #this means the segment in in the lower arm of the chromosome
					q_arm_len += segmentLength
						
					if int(outputArray[10]) > 0 : #cnt is more than 0 (if cnt is 0, this is either a deletion or no coverage, either way, it is not indicative of LOH event)
						if outputArray[11] == "0": #estimated number of A alleles is 0
							q_arm_lossA = segmentLength
							chr_lossA += segmentLength
						if outputArray[12] == "0": #estimated number of B alleles is 0
							q_arm_lossB += segmentLength
							chr_lossB += segmentLength


		####CALCULATE LOSSES FOR EACH ARM OF CHROM (..)#######
		# We discard ANY chromosome or ANY chromosomal arm that has 90% or more loss as these do not indicate LOH results
		#Case 1: both arms have less than 90% loss and get appended to valid chromosomes
		#Case 2: both arms have more than 90% loss and do not get appended to valid chromosomes
		#Case 3: only one of the arms have a loss of less than 90% and gets appened to valid chromosomes

		if not chr_len == 0 :
			percentAloss_chr = float(chr_lossA) / chr_len
			percentBloss_chr = float(chr_lossB) / chr_len
			print "====================================================================="
			print "Aloss over entire chrom " + str(chromosome) + ": " + str(percentAloss_chr)
			print "Bloss over entire chrom " + str(chromosome) + ": " + str(percentBloss_chr)

			#Check that entire loss over chromosome for allela A/allele B < .90
			if percentAloss_chr < float(loss_threshold) and percentBloss_chr < float(loss_threshold):

				if not p_arm_len == 0:
					percentAloss_p_arm = float(p_arm_lossA) / p_arm_len
					percentBloss_p_arm = float(p_arm_lossB) / p_arm_len
					
					#if so, check that that loss over the p arm for allela A/allele B < .90
					if percentAloss_p_arm < float(loss_threshold) and percentBloss_p_arm < float(loss_threshold): #CHECKKKKKKKK if this is and
						validChromosomeArms = validChromosomeArms + str(chromosome) + "_p_arm,"
						total_bases  += p_arm_len
						print "p arm loss for allele A: " + str(p_arm_lossA)
						print "p arm loss for allele B: " + str(p_arm_lossB )
						total_lost_bases += p_arm_lossA + p_arm_lossB
						print "total lost bases after chr "+ str(chromosome) + " arm p is " + str(total_lost_bases)
				
				if not q_arm_len == 0:
					percentAloss_q_arm = float(q_arm_lossA) / q_arm_len
					percentBloss_q_arm = float(q_arm_lossB) / q_arm_len

					#and, check that that loss over the q arm for allela A/allele B < .90 
					if percentAloss_q_arm < float(loss_threshold) and percentBloss_q_arm < float(loss_threshold):
						validChromosomeArms = validChromosomeArms + str(chromosome) + "_q_arm,"
						total_bases  += q_arm_len
						print "q arm loss for allele A: " + str(q_arm_lossA)
						print "q arm loss for allele B: " + str(q_arm_lossB)
						total_lost_bases += q_arm_lossA + q_arm_lossB
						print "total lost bases after chr "+ str(chromosome) + " arm q is " + str(total_lost_bases)
				
					

def write_to_output_file(output_file):
	with open(output_file, "w+") as output:
		output.write("Total valid bases analyzed: %d\n" %total_bases)
		output.write("Valid lost bases identified: %d\n" %total_lost_bases)
		lossScore = 0

		if not total_lost_bases == 0:
			lossScore = float(total_lost_bases)/total_bases

		print "loss score is " + str(lossScore)
		output.write("Fraction valid bases lost: %.10f\n" %lossScore)
		output.write("Valid chromosomal arms: %s" %str(validChromosomeArms).strip(","))



#MAIN CODE#
centromere_file =  sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]
loss_threshold = sys.argv[4]


def LOH_analysis (centromere_file, input_file, output_file, loss_threshold):
	#Find out centromere info
	centromereInfo(centromere_file)

	#check that copyno file is modified if a segment in it overlaps with a centromere
	modified_input_file = input_file.replace(".txt", "_modified.txt") 
	checkCopyNoFile (input_file, modified_input_file)
	
	#Calculate LOH score
	for chromosome in chromosomesToTest:
		calculateLOHscore(modified_input_file, chromosome, loss_threshold)

	#write results to file
	write_to_output_file(output_file)


LOH_analysis (centromere_file, input_file, output_file, loss_threshold)

#For my own reference
print "\n"
print "Valid Chromosomal Arms: " + str(validChromosomeArms).strip(",")

print "total_bases: " + str(total_bases)
print "total_lost_bases: " + str(total_lost_bases)
