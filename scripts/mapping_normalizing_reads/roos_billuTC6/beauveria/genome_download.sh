#!/bin/bash
#SBATCH -J genome_downloads_beau
#SBATCH -t 4:0:0
#SBATCH --mem=5G
#SBATCH -o /home/uu_bio_fg/rbrouns/data/sbatch_out/TC6_genome_downloads_out.out
#SBATCH -e /home/uu_bio_fg/rbrouns/data/sbatch_out/TC6_genome_downloads_out.e
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=r.c.a.brouns2@students.uu.nl

### Download the genome and annotation file from NCBI
## The accession number is  GCA_000280675.1_ASM28067v1
## IDs: 411938 [UID] 411938 [GenBank] 1202398 [RefSeq]

# 00: Set directory
cd /home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/genome/beau/

# 01: Download the genome file 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/280/675/GCA_000280675.1_ASM28067v1/GCA_000280675.1_ASM28067v1_genomic.fna.gz
# unzip the genome file
gzip -d GCA_000280675.1_ASM28067v1_genomic.fna.gz
# rename the genome file 
mv GCA_000280675.1_ASM28067v1_genomic.fna beau_ARSEF2860_genome.fna

# 02: Download the annotation file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/280/675/GCA_000280675.1_ASM28067v1/GCA_000280675.1_ASM28067v1_genomic.gff.gz
# unzip the annotation file
gzip -d GCA_000280675.1_ASM28067v1_genomic.gff.gz
# rename the annotation file
mv GCA_000280675.1_ASM28067v1_genomic.gff beau_ARSEF2860_genome.gff

### END