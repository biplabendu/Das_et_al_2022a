#!/bin/bash
#SBATCH --cpus-per-task=8


#load modules (packages)
#module load bbmap # for bbduk.sh

#set directory (raw seq reads)
cd /home/billu/TC6/data/deBekker2017_ophio/raw_reads/GSE101312/


### STEP: Trimming 
#
# TO DO: ADD FOR-LOOP TO INCLUDE ALL SAMPLES?
#
## Adapter trim and quality trim in one run
#bbduk.sh -Xmx1g in1=TC6-02A_ATCACG.fastq.gz out1=./trimmed_reads/trimmed_TC6-02A_ATCACG.fastq.gz minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 literal=ATCACG,CGATGT hdist=1 tpe tbo
#
# "gtrim=10" will quality trim at Q10 using Phred algorithm, "qtrim=rl" will trim right and left sides
# "ref=truseq.fq.gz" included Illumina truseq adapters within BBMap package
# "hdist=1" allows one mismatch, "k=25" is k-mer size of 25, "mink=11" allow shorter k-mers of 11 at end of read
# default looks at reverse-complement and forward seuence of refseq otherwise "rcomp=<t/r>"
# "-Xmx1g" flag tells BBDuk to use 1GB of RAM
#
#
## Or split into adapter trimming first
#
#for file in *fastq
#do
#
#bbduk.sh \
#	-Xmx1g \
#	in1=${file} \
#	out1=./trimmed_reads/trimmed_A_${file} \
#	ref=/share/apps/bbmap/resources/adapters.fa \
#	#literal=ATCACG,CGATGT \
#	ktrim=r \
#	k=23 \
#	mink=11 \
#	hdist=1 \
#	tpe \
#	tbo	
# "tbo" (trim adapters based on pair overlap detection) and "tpe" (both reads to the same length) flags for normal paired-end fragment libraries
#
## Then quality trimming
#
#bbduk.sh \
#	-Xmx1g \
#	in=./trimmed_reads/trimmed_A_${file} \
#	out=./trimmed_reads/trimmed_AQ_${file} \
#	qtrim=rl \
#	trimq=10
#done
#
# remove all the trimmed_A_* files
#rm ./trimmed_reads/trimmed_A_*
#
### chunk ends


### STEP: Run FASTQC on the trimmed reads
#
# specify directory of the *fastq.gz files
cd ./trimmed_reads/
#
# run fastQC on all files and save output in the "output_fastqc" folder
fastqc -o ./output_fastqc *fastq
#
#
### chunk ends
