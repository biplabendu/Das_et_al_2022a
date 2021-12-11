#!/bin/bash
#SBATCH --cpus-per-task=8


## PART ONE: Get the protein fasta files

# 00: Set directory
cd /home/billu/TC6/data/protein_fasta

###########
## OPHIO ##
###########
# 01: Download the protein fasta file 
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_protein.faa.gz
# unzip the file
#gzip -d GCA_012980515.1_Ophcf2_protein.faa.gz
# rename the file 
#mv GCA_012980515.1_Ophcf2_protein.faa ophcf2_protein

##########
## BEAU ##
##########
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/280/675/GCF_000280675.1_ASM28067v1/GCF_000280675.1_ASM28067v1_protein.faa.gz
#gzip -d GCF_000280675.1_ASM28067v1_protein.faa.gz
#mv GCF_000280675.1_ASM28067v1_protein.faa beau_protein
#
## DONE.
#
## part one ends
# ----------------------------------------------------------------


## PART TWO: Run proteinortho on the protein fasta files

# o.cflo vs b.bassiana 
#################
proteinortho5 -step=1 -project=ophcf2_v_beau_1 -p=blastp+ ./beau_protein ./ophcf2_protein
## done
proteinortho5 -step=2 -project=ophcf2_v_beau_1 -p=blastp+ ./beau_protein ./ophcf2_protein
#
##
#
## part two ends
# ----------------------------------------------------------------


