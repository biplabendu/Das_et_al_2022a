#!/bin/bash
#SBATCH --cpus-per-task=8

## set directory
# for Roos
#cd /home/fg/r_brouns/ant_fungus/TC6/data/raw_seq_reads/
#
# for Billu
#cd /home/billu/TC6/data/raw_seq_reads/

## set path for the data
# for Roos
#path = /home/fg/r_brouns/ant_fungus/TC6/data/
# for Billu
path=/home/billu/TC6/data/

### STEP: Trimming 
#
# input file directory: /home/billu/TC6/data/raw_seq_reads/
# input file names: TC6-02A_ATCACG.fastq.gz (example)
#
# [adapter trimming + quality trimming]
#
# output file names: trimmed_AQ_TC6-02A_ATCACG.fastq.gz (example)
#
### chunk ends


### 01_Code_Chunk: Aligning reads 
#
# TO DO: ALTER FOR FUNGAL GENOME, BOTH OPHIO AND BBAS.
#
# create a directory for storing the genome and annotation files
#cd $path/genome/ophio/
#
# Things to download from NCBI webstite, unzip, and rename:
#
# 01: Genome: Ophcf2_genome.fna
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.fna.gz
#gzip -d GCA_012980515.1_Ophcf2_genomic.fna.gz
#mv GCA_012980515.1_Ophcf2_genomic.fna ophcf2_genome.fna
#
# 02: Annotation file: Ophcf2_genome.gff
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.gff.gz
#gzip -d GCA_012980515.1_Ophcf2_genomic.gff.gz
#mv GCA_012980515.1_Ophcf2_genomic.gff ophcf2_genome.gff
#
## index the Ophio genome, in parent dir (seems redundant)
# hisat2-build -f ophcf2_genome.fna ophcf2_index 
# "-f" for input are FASTA files
#
# fresh ophio index, with exons
#gffread "ophcf2_genome.gff" -T -o ophcf2.gtf
#hisat2_extract_splice_sites.py "ophcf2.gtf" > ophcf2_splicesites.txt
#hisat2_extract_exons.py "ophc2.gtf" > ophcf2_exons.txt
#
# echo "BUILDING THE OPHIO INDEX WITH SPLICE AND EXONS SITE"
# index the exon site version of the ophcf2 genome
#hisat2-build -f --ss ophcf2_splicesites.txt --exon ophcf2_exons.txt ophcf2_genome.fna ophcf2_exons_index
# "-ss" is for splice sites, followed by "--exon" exon sites
#
# map (trimmed) rnaseq sample to (indexed) genome
#
### TO DO: MAKE A LOOP TO INDIVIDUALLY MAP ALL RNASEQ SAMPLES TO THE GENOME
#
cd $path/trimmed_reads/ophio/
for file in *fastq.gz
do
hisat2 \
	-p 8 
	--new-summary \
	-q \
	--dta-cufflinks \
	-x $path/genome/ophio/ophcf2_exons_index \
	-U $path/trimmed_reads/ophio/${file} \
	-S $path/mapped_reads/ophio/mapped_${file%%.*}.sam
done
#
### chunk ends

