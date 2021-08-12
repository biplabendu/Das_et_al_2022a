#!/bin/bash
#SBATCH --cpus-per-task=8


#load modules (packages)
module load hisat2
module load cufflinks
module load samtools
module load bbmap

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
path_data=/home/billu/TC6/data/deBekker2017_ophio/

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
cd $path_data/genome/
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
# hisat2-build -f GCA_001272575.2_ASM127257v2_genomic.fna ophio_kim_index
# "-f" for input are FASTA files
#
# fresh ophio index, with exons
#gffread "GCA_001272575.2_ASM127257v2_genomic.gff" -T -o ophio_kim.gtf
#hisat2_extract_splice_sites.py "ophio_kim.gtf" > ophio_kim_splicesites.txt
#hisat2_extract_exons.py "ophio_kim.gtf" > ophio_kim_exons.txt
#
#echo "BUILDING THE OPHIO INDEX WITH SPLICE AND EXONS SITE"
# index the exon site version of the ophcf2 genome
#hisat2-build \
#	-f \
#	--ss ophio_kim_splicesites.txt \
#	--exon ophio_kim_exons.txt \
#	GCA_001272575.2_ASM127257v2_genomic.fna \
#	ophio_kim_exons_index
#
# map (trimmed) rnaseq sample to (indexed) genome
#
### TO DO: MAKE A LOOP TO INDIVIDUALLY MAP ALL RNASEQ SAMPLES TO THE GENOME
#
cd $path_data/raw_reads/GSE101312/trimmed_reads/
for file in *fastq
do
hisat2 \
	-p 8 \
	--new-summary \
	-q \
	--dta-cufflinks \
	-x $path_data/genome/ophio_kim_exons_index \
	-U $path_data/raw_reads/GSE101312/trimmed_reads/${file} \
	-S $path_data/mapped_reads/mapped_${file%%.*}.sam
done
#
### chunk ends

