# HRD Testing Project

Important links:
 - Tracker issue in genetics/projects repo: https://gitlab.labmed.uw.edu/genetics/projects/issues/100
 - Current research pipeline: https://bitbucket.org/nithishak/loh-pipeline/src/master/

## Creating hg19.centromere.txt
This file is dependent on genome build (e.g.: hg19 vs hg38)
    - go to the genome browser: https://genome.ucsc.edu/cgi-bin/hgTables
    - select relevant genome build from "assembly" dropdown
    - select `all tables` from "group" dropdown
    - select `gap` from "table" dropdown
    - select `genome` from "region" bullets
    - click `get output`
This file gives the coordinates of genomic features (including centromere).
The LOH scoring script parses this file to decide which segments are in the P vs Q chromosome arms.

## Pipeline installation

Install the nextflow binary in this directory
  
```bash
wget -qO- https://get.nextflow.io | bash
```
  
Execute locally, using docker images (must be available via docker locally)
  
```bash
./nextflow run main.nf --tumor tumor_dir/ --normal normal_dir/ --gc hg19.wig --cen hg19.centromere.txt -profile docker
```
