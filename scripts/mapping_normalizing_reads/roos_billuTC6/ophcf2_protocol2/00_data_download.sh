#!/bin/bash
#SBATCH -J data_download 
#SBATCH -t 4:0:0
#SBATCH --mem=5G


### Download Data
# This scripts downloads the raw RNA-seq reads, along with the genome and annotation files of Ophio camp-florani. 

## 00: Give the path to directory of the experiment_name
path = /home/set/absolute/path/to/directory #don’t add a / at the end of the path
experiment_name = TC6



## TO DO: write script for download data for RNA-seq data
# The raw_read are in .fastq.gz format


## Code to download the genome and annotation file from NCBI website. The genome and annotation file are necessary to map the reads. The genome file should be .fna format and the annotation file .gff format.
#
# 01: Make directory to store genome and annotation file and set to work from that directory
mkdir $path/$experiment_name/data/genome
cd $path/$experiment_name/data/genome/
#
# 02: Download the genome file: Ophcf2_genome.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.fna.gz
# unzip the genome file
gzip -d GCA_012980515.1_Ophcf2_genomic.fna.gz
# rename the genome file 
mv GCA_012980515.1_Ophcf2_genomic.fna ophcf2_genome.fna
#
# 03: Download the annotation file: Ophcf2_genome.gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.gff.gz
# unzip the annotation file
gzip -d GCA_012980515.1_Ophcf2_genomic.gff.gz
# rename the annotation file
mv GCA_012980515.1_Ophcf2_genomic.gff ophcf2_genome.gff
#
#
### END
