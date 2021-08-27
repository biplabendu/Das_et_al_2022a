#!/bin/bash
#SBATCH -J trim_reads
#SBATCH -t 4:0:0
#SBATCH --mem=5G

### Trim reads
# This scripts trims the raw reads for adapetrs and quality with the bbduk tool. The trimmed reads are stored in a trimmed_reads output directory. 
# After the reads are trimmed, it will run fastqc, to check the quality of the reads. This output is stored in a fastqc_output directory.
# This scripts needs the software of bbduk and fastqc, and make sure that it is working
# The raw reads are stored in the folder /home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/ophio/demultiplexed_reads/ 

# Code to trim the reads
#
# 00: Make output folder for the trimmed reads
path = /home/uu_bio_fg/rbrouns/data/ant_fungus # do not end with /
experiment_name = TC6
mkdir $path/$experiment_name/data/ophio/trimmed_reads
#
# set directory (for raw seq reads)
cd /$path/$experiment_name/data/ophio/demultiplexed_reads/

### Trimming 
#
## 01: First adapter trimming of all reads
# For this bit of code the file with the used adapeters is needed at the parameter ref=<file.fa>. 
# The path used in this script is al follows, /home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/bbmap_adapters/adapters.fa
#
for file in *fastq.gz
do
#
bbduk.sh \
	-Xmx1g \
	in1=${file} \
	out1=./../ophio/trimmed_reads/trimmed_A_${file} \
	ref=/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/bbmap_adapters/adapters.fa \
	#literal=ATCACG,CGATGT \
	ktrim=r \
	k=23 \
	mink=11 \
	hdist=1 \
	tpe \
	tbo	

# Parameter explanation:
	#"gtrim=10" will quality trim at Q10 using Phred algorithm, "qtrim=rl" will trim right and left sides
	# "ref=truseq.fq.gz" included Illumina truseq adapters within BBMap package
	# "hdist=1" allows one mismatch, "k=25" is k-mer size of 25, "mink=11" allow shorter k-mers of 11 at end of read
	# default looks at reverse-complement and forward seuence of refseq otherwise "rcomp=<t/r>"
	# "-Xmx1g" flag tells BBDuk to use 1GB of RAM
	# "tbo" (trim adapters based on pair overlap detection) and "tpe" (both reads to the same length) flags for normal paired-end fragment libraries
#
## 02: Then quality trimming of all reads that already have been trimmed for adapaters
#
bbduk.sh \
	-Xmx1g \
	in=./../ophio/trimmed_reads/trimmed_A_${file} \
	out=./../ophio/trimmed_reads/trimmed_AQ_${file} \
	qtrim=rl \
	trimq=10
done
#
# remove all the trimmed_A_* files, which are the files that are only trimmed on adapters. What remains are the adapter and quality trimmed reads.
rm ./../trimmed_reads/ophio/trimmed_A_*
#
### END


### 04: Run FastQ on the trimmed reads
#
# Specify directory of the *fastq.gz files, which are the trimmed reads
cd /$path/$experiment_name/data/ophio/trimmed_reads/
# Create a output_fastqc folder 
mkdir fastqc_output
#
# Run fastQC on all files and save output in the "fastqc_output" folder
fastqc -o ./fastqc_output/ *fastq.gz
#
### END
