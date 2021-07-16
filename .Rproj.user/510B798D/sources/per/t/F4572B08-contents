#!/bin/bash
#SBATCH --cpus-per-task=8


#load modules (packages)
module load hisat2
module load cufflinks
module load samtools


### STEP: Trimming 
#
# TO DO: ADD A CODE CHUNK HERE TO ELIMINATE ADAPTER SEQUENCES AND LOW-QUALITY READS
# (use Trimmomatic or Cutadapt or BBDuk)
#
### chunk ends


### 01_Code_Chunk: Aligning reads 
#
# index the ant genome, in parent dir
hisat2-build -f camp_genome.fna camp_index
#
# fresh ant index, with exons
gffread "/home/billu/cflo_genome/GCF_003227725.1_Cflo_v7.5_genomic.gff" -T -o camp.gtf
hisat2_extract_splice_sites.py "/home/billu/cflo_genome/camp.gtf" > camp_splicesites.txt
hisat2_extract_exons.py "/home/billu/cflo_genome/camp.gtf" > camp_exons.txt
#
# echo "BUILDING THE ANT INDEX WITH SPLICE AND EXONS SITE"
# index the exon site version of the ANT genome
hisat2-build -f --ss camp_splicesites.txt --exon camp_exons.txt ./camp_genome.fna camp_exons_index
#
# map timecourse rnaseq samples to (indexed) genome
#
### TO DO: MAKE A LOOP TO INDIVIDUALLY MAP ALL RNASEQ SAMPLES TO THE GENOME
#
hisat2 -p 8 --new-summary -q --dta-cufflinks -x /home/billu/cflo_genome/camp_exons_index -U /home/billu/TC5_seq_reads/Foragers/Trimmed-2F_Galaxy96.fastq.gz -S /home/billu/TC5_seq_reads/2F_time5_1.sam
#
### chunk ends

### 02_Code_Chunk: Counting and Normalization (MIGHT NEED TO CHANGE)
#
# using .bam from hisat2 mapping>samtools sort and convert 
# NO CUFFLINKS step used (nor cuffquant).
# -b ant genome
# sort and convert Hisat2 sams to sorted bams
samtools sort -o ./2F_time5_1.bam ./2F_time5_1.sam
samtools sort -o ./4F_time5_1.bam ./4F_time5_1.sam
#
# use the following code to obtain normalized FPKM values
# TO DO: Make a loop to obtain normalized gene counts for each sample
#
cuffdiff -o Foragers -b "/home/billu/cflo_genome/camp_genome.fna" -p 16 -u "/home/billu/cflo_genome/GCF_003227725.1_Cflo_v7.5_genomic.gff" \2F_time5_1.bam \4F_time5_1.bam
#
# TO DO: Check HTSeq-based TMM normalization instead of Cuffdiff-based FPKM normalization
#
### chunk ends