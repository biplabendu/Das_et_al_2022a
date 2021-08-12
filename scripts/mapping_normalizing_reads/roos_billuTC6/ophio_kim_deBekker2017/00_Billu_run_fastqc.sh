#!/bin/bash
#SBATCH --cpus-per-task=8

# specify directory of the *fastq.gz files
cd /home/billu/TC6/data/deBekker2017_ophio/raw_reads/GSE101312/

# run fastQC on all files and save output in the "output_fastqc" folder
fastqc -o ./analysis/output_fastqc *fastq
