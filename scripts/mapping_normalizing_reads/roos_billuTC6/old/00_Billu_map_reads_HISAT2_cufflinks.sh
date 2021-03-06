#!/bin/bash
#SBATCH --cpus-per-task=8


#load modules (packages)
module load hisat2
module load cufflinks
module load samtools
module load bbmap

# set directory
# for Roos
cd /home/fg/r_brouns/ant_fungus/TC6/data/raw_seq_reads/
# for Billu
# cd /home/billu/TC6/data/raw_seq_reads/

# set path
path = /home/fg/r_brouns/ant_fungus/TC6/data/

### STEP: Trimming 
#
# TO DO: AD LOOP FOR TO INCLUDE ALL SAMPLES?
#
## Adapter trim and quality trim in one run
bbduk.sh -Xmx1g in1=<r1.fq> in2=<r2.fq> out1=<clean1.fq> out2=<clean2.fq> minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 ref=truseq.fq.gz hdist=1
# "gtrim=10" will quality trim at Q10 using Phred algorithm, "qtrim=rl" will trim right and left sides
# "ref=truseq.fq.gz" included Illumina truseq adapters within BBMap package
# "hdist=1" allows one mismatch, "k=25" is k-mer size of 25, "mink=11" allow shorter k-mers of 11 at end of read
# default looks at reverse-complement and forward seuence of refseq otherwise "rcomp=<t/r>"
# "-Xmx1g" flag tells BBDuk to use 1GB of RAM
#
## Or split into adapter trimming first
bbduk.sh -Xmx1g in=<reads.fq> out=<clean_A.fq> ref=<adapters.fa> ktrim=r k=23 mink=11 hdist=1 tpe tbo
# "tbo" (trim adapters based on pair overlap detection) and "tpe" (both reads to the same length) flags for normal paired-end fragment libraries
#
## Then quality trimming
bbduk.sh -Xmx1g in=<clean_A.fq> out=<clean_AQ.fq> qtrim=rl trimq=10
#
### chunk ends


### 01_Code_Chunk: Aligning reads 
#
# TO DO: ALTER FOR FUNGAL GENOME, BOTH OPHIO AND BBAS.
#
# Things to download from NCBI webstite, unzip, and rename:
# 01: Genome: Ophcf2_genome.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.fna.gz
gzip -d GCA_012980515.1_Ophcf2_genomic.fna.gz
mv GCA_012980515.1_Ophcf2_genomic.fna ophcf2_genome.fna
# 02: Annotation file: Ophcf2_genome.gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.gff.gz
gzip -d GCA_012980515.1_Ophcf2_genomic.gff.gz
mv GCA_012980515.1_Ophcf2_genomic.gff ophcf2_genome.gff
#
## index the Ophio genome, in parent dir (seems redundant)
# hisat2-build -f ophcf2_genome.fna ophcf2_index 
# "-f" for input are FASTA files
#
# fresh ophio index, with exons
gffread "ophcf2_genome.gff" -T -o ophcf2.gtf
hisat2_extract_splice_sites.py "ophcf2.gtf" > ophcf2_splicesites.txt
hisat2_extract_exons.py "ophc2.gtf" > ophcf2_exons.txt
#
# echo "BUILDING THE OPHIO INDEX WITH SPLICE AND EXONS SITE"
# index the exon site version of the ANT genome
hisat2-build -f --ss ophcf2_splicesites.txt --exon ophcf2_exons.txt ophcf2_genome.fna ophcf2_exons_index
#
# map timecourse rnaseq samples to (indexed) genome
#
### TO DO: MAKE A LOOP TO INDIVIDUALLY MAP ALL RNASEQ SAMPLES TO THE GENOME
#
hisat2 -p 8 --new-summary -q --dta-cufflinks -x ophcf2_exons_index -U ${path}/trimmed_reads/trimmed_AQ_TC6-02A_ATCACG.fastq.gz -S ${path}/mapped_TC6-02A_ATCACG.sam
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
# TO DO: Check HTSeq-based TMM normalization instead of Cuffdiff-based FPKM/RPKM normalization
#
### chunk ends
