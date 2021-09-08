#!/bin/bash
#SBATCH -J mapping
#SBATCH -t 4:0:0
#SBATCH --mem=5G

### Mapping 
# This script maps the reads to the reference genome, with HISAT2. 
# For this script the software of HISAT2 and gffread is needed.
#
# 01: Index the O. camp-florani genome, with exon and and splice sites. 
#
# Set genome folder as working directory to store the index files we're creating
path = /home/uu_bio_fg/rbrouns/data/ant_fungus # don't ad / at the end
experiment_name = TC6
cd $path/$experiment_name/data/ophio/genome
#
# Since we work with RNA-seq data (and thus only transcibed DNA, which are introns), we want to index the exons and slice sites from the mapping.
gffread "ophcf2_genome.gff" -T -o ophcf2.gtf 
hisat2_extract_splice_sites.py "ophcf2.gtf" > ophcf2_splicesites.txt
hisat2_extract_exons.py "ophcf2.gtf" > ophcf2_exons.txt
#
# index the exon site version of the ophcf2 genome with hisat2-build
hisat2-build \ 
    -f \ # The reference input files are FASTA files
    --ss ophcf2_splicesites.txt \ # Input the splice sites
    --exon ophcf2_exons.txt \ # Input the exons
    ophcf2_genome.fna \ # Input the reference genome
    ophcf2_exons_index # Basename of output file to write indexed genome with exons and splice sites in


# 02: Map the trimmed RNAseq sample reads to indexed genome
#
# Make a folder to store the mapped reads
mkdir $path/$experiment_name/data/ophio/mapped_reads
#
# Set folder with trimmed reads as working directory
cd $path/$experiment_name/data/ophio/trimmed_reads/
#
# Map the reads with hisat2 by looping over each file in the folder
for file in *fastq.gz
do
hisat2 \
        -p 8 \ # use 8 cpu's
        --new-summary \ # Print alignment summary in a new style, which is more machine-friendly.
        -q \ # Reads are FASTQ files.
        --dta-cufflinks \ # Report alignments tailored for transcript assembler Cufflinks
        -x $path/$experiment_name/data/ophio/genome/ophcf2_exons_index \ # Input basename of index genome file
        -U $path/$experiment_name/data/ophio/trimmed_reads/${file} \ # Input the file containing unpaired reads to be aligned (can also make a comma seperated list of the files)
        -S $path/$experiment_name/data/ophio/mapped_reads/mapped_${file%%.*}.sam # File to write SAM alignments to
done
#
### END