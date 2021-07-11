#!/bin/bash -l
# author: rfitak
#SBATCH -J billu
#SBATCH -o billu.out
#SBATCH -e billu.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 2:00:00
#SBATCH --mem-per-cpu=16000

# NOTES:
# -J  :: name of job
# -o  :: name of file to write job output (STDOUT)
# -e  :: name of file to write job error  (STDERR)
# -p  :: partition or queue to run job, no need to change this on coombs
# --cpus-per-task=1   :: how many cpus do you need?  Max is 32
# --mem-per-cpu=8000  :: how much memory is needed for each cpu (in MB)? Max = 16000

# Begin normal BASH commands
cd /home/billu/TC6/data/deBekker2017_ophio/genome/  # (BILLU: change this to your folder)

# Download Ophiocordyceps kimflemingiae genome and annotation files from NCBI
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/272/575/GCA_001272575.2_ASM127257v2/GCA_001272575.2_ASM127257v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/272/575/GCA_001272575.2_ASM127257v2/GCA_001272575.2_ASM127257v2_genomic.gff.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.gff.gz
gunzip GCA_001272575.2_ASM127257v2_genomic.fna.gz
gunzip GCA_001272575.2_ASM127257v2_genomic.gff.gz

# TO DO: