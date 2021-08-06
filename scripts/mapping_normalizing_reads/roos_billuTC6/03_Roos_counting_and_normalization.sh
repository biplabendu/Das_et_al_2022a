#!/bin/bash
#SBATCH -J counting_normalization_ophio 
#SBATCH -t 4:0:0
#SBATCH --mem=5G
#SBATCH -o /home/uu_bio_fg/rbrouns/data/sbatch_out/TC6_ophio_couting_normalization_out.out
#SBATCH -e /home/uu_bio_fg/rbrouns/data/sbatch_out/TC6_ophio_counting_normalization_out.e
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=r.c.a.brouns2@students.uu.nl


### STEP: Counting and Normalization
#
# define path to mapped reads
path=/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/

# using .bam from hisat2 mapping>samtools sort and convert 
# NO CUFFLINKS step used (nor cuffquant).
# sort and convert Hisat2 sams to sorted bams
for file in $path/mapped_reads/ophio/*.sam
do
samtools sort -o ${file%%.*}.bam ${file}
done

### Cuffdiff normalization and counting
#
# Make directoty for output cuffdiff
mkdir $path/cuffdiff_out
mkdir $path/cuffdiff_out/ophio

## Option 1: Single sample vs sample for cuffdiff
cd ${path}/mapped_reads/ophio/
#
sample1=mapped_trimmed_AQ_TC6-02A_ATCACG.bam
sample2=mapped_trimmed_AQ_TC6-02A_ATCACG.bam
#
# Run cuffdiff 
cuffdiff \
        -o ${path}/cuffdiff_out/ophio/ \
        -b ${path}/genome/ophio/ophcf2_genome.fna \
        -p 16 \
        -u ${path}/genome/ophio/ophcf2_genome.gff \
        ${sample1} \
        ${sample2}


## Option 2: Nested loop for pairwise comparison of TC2 vs TC4, TC2 vs TC6 ect., TC4 vs TC6, TC4 vs TC8, ect.
# Set directory
cd $path/mapped_reads/ophio/
#
for sample1 in *.bam
do
    for sample2 in *.bam
    do
        cuffdiff \
            -o $path/cuffdiff_out/ophio \
            -b "/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/genome/ophio/ophcf2_genome.fna" \
            -p 16 \
            -u "/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/genome/ophio/ophcf2_genome.gff" \
            $sample1 \
            $sample2
    done
done


## Option 3: Loop of only sample 1 vs all (TC2 vs all)
cd ${path}/mapped_reads/ophio/
sample1=mapped_trimmed_AQ_TC6-02A_ATCACG.bam
#
# Run cuffdiff for sample1 vs all
for sample2 in *.bam
do
        cuffdiff \
                -o ${path}/cuffdiff_out/ophio/ \
                -b ${path}/genome/ophio/ophcf2_genome.fna \
                -p 16 \
                -u ${path}/genome/ophio/ophcf2_genome.gff \
                ${sample1} \
                ${sample2}
done


### Rename file with gene expression for sample vs sample
#
# first make variable for pairwise name of file
pw_name1=$(echo "$sample1" | cut -d '-' -f 2 | cut -d '_' -f 1)
pw_name2=$(echo "$sample2" | cut -d '-' -f 2 | cut -d '_' -f 1)
# rename file
mv ${path}/cuffdiff_out/ophio/gene_exp.diff ${path}/cuffdiff_out/ophio/gene_exp_${pw_name1}vs${pw_name2}.csv
# remove all the other files 
# (does not work on cluster of UU_
cd ${path}/cuffdiff_out/ophio/
rm !(*.csv)   

## Only put cols with gene expression per sample in a csv
#
# first make a temporary csv file with only the cols with test_id, gene_id, gene, and locus
csv_name=TC6_ophio
cut -d$'\t' -f 1,2,3,4 gene_exp_${pw_name1}vs${pw_name1}.csv > ${csv_name}_temp_start.csv
# Then paste 8 column from pairwise comparison on sample into new csv file
paste ${csv_name}_temp_start.csv <(awk '{print $8}' gene_exp_${pw_name1}vs${pw_name2}) > ${csv_name}_gene_exp.csv
# Rename header of column with value of sample from value_2 to sample name
header_sample=TP_${pw_name2}


# append with -f 8 from the other csv files
# work in progress

#
# TO DO: Check HTSeq-based TMM normalization instead of Cuffdiff-based FPKM/RPKM normalization
#
### chunk ends