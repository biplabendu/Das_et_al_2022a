#!/bin/bash -l
# author: billu
#SBATCH -J ophio_kim_get_data
#SBATCH -o ophio_kim_get_data.out
#SBATCH -e ophio_kim_get_data.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=biplabendu.das@knights.ucf.edu

# NOTES:
# -J  :: name of job
# -o  :: name of file to write job output (STDOUT)
# -e  :: name of file to write job error  (STDERR)
# -p  :: partition or queue to run job, no need to change this on coombs
# --cpus-per-task=1   :: how many cpus do you need?  Max is 32
# --mem-per-cpu=8000  :: how much memory is needed for each cpu (in MB)? Max = 16000


# CODE CHUNK BEGINS
#########################################
## DOWNLOAD GENOME AND ANNOTATION FILE
#########################################
#
# Set your working directory
#cd /home/billu/TC6/data/deBekker2017_ophio/genome/  # (BILLU: change this to your folder)
#
# Download genome from NCBI
#
# Ophiocordyceps-camponoti-floridani
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.fna.gz
#
# Ophiocordyceps-kimflemingiae
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/272/575/GCA_001272575.2_ASM127257v2/GCA_001272575.2_ASM127257v2_genomic.fna.gz
#
# Download the annotation file from NCBI
# 
# Ophiocordyceps-camponoti-floridani
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.gff.gz
#
# Ophiocordyceps-kimflemingiae
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/272/575/GCA_001272575.2_ASM127257v2/GCA_001272575.2_ASM127257v2_genomic.gff.gz
#
#
# Unzip the genome and annotation files
#gunzip GCA_001272575.2_ASM127257v2_genomic.fna.gz
#gunzip GCA_001272575.2_ASM127257v2_genomic.gff.gz
#
#########################################
# CODE CHUNK ENDS


# CODE CHUNK BEGINS
#########################################
## DOWNLOAD RAW SEQUENCING FILES
#########################################
#
# Set your working directory
cd /home/billu/TC6/data/deBekker2017_ophio/raw_reads/GSE101312/  
#
# Download all raw sequencing files for accession number GSE101312 (de Bekker 2017)
#
# load the sratoolkit
module load sratoolkit
#
for i in $(cat SRR_Acc_List_GSE101312.txt); do fastq-dump $i; done
#

