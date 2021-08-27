#!/bin/bash
#SBATCH -J counting_normalization
#SBATCH -t 4:0:0
#SBATCH --mem=5G

### Counting and Normalization
# This script count the mapped reads and normimlizes them to get the gene expression values in FRKPM
# This script requires Samtools and cuffdiff (from CuffLinks)

# 00: define path to directory
experiment_name = TC6
path = /home/uu_bio_fg/rbrouns/data/ant_fungus/$experiment_name

## 01: Sort and convert all Hisat2 sams files to sorted bams with samtools sort
# 
cd $path/data/ophio/mapped_reads/
#
for file in *.sam
do
samtools sort \
    -o ${file%%.*}.bam \ # Output same basename with .bam extension
    ${file} # Input .sam files
done

## 02: Cuffdiff normalization
# Cuffdiff compare one sample to another to obtain normalized gene expression values. Therefore we define our first sample and compare all the other samples to this sample. 
# Per sample vs sample analysis we create a output file with the gene expressoin values.
#
# Create output folder
mkdir $path/data/ophio/cuffdiff_out
#
# Set working directory to mapped reads
cd ${path}/data/ophio/mapped_reads/
#
# 
# Name the first sample (which has 2 in the file name because first sample is timepoint 2)
sample1=mapped_trimmed_AQ_*02*.bam
#
# Run cuffdiff for sample1 vs all samples
for sample2 in *.bam
do
        cuffdiff \
                -o ${path}/data/ophio/cuffdiff_out/ \ # Output folder
                -b ${path}/data/ophio/genome/ophcf2_genome.fna \ # Input reference genome
                -p 16 \ # Number of threads
                -u ${path}/data/ophio/genome/ophcf2_genome.gff \ # Input annotation file
                ${sample1} \ # run cuffdiff on sample 1
                ${sample2} # vs sample 2

        # Then rename file with gene expression to get per analysis a output file
        # make variable for pairwise name of file
        pw_name1=$(echo $sample1 | cut -d '-' -f 2 | cut -d "_" -f 1)
        pw_name2=$(echo "$sample2" | cut -d '-' -f 2 | cut -d '_' -f 1)
        # rename file to excuted analysis
        mv ${path}/data/ophio/cuffdiff_out/gene_exp.diff ${path}/data/ophio/cuffdiff_out/gene_exp_${pw_name1}vs${pw_name2}.csv
done

# remove all the other files 
cd ${path}/data/ophio/cuffdiff_out/
rm !(*.csv) 
#
### END  