#!/bin/bash
#SBATCH -J counting_normalization_deBekker2017
#SBATCH -t 4:0:0
#SBATCH --mem=5G

module load cufflinks
module load samtools

### STEP: Counting and Normalization
#
# define path to mapped reads
path_data=/home/billu/TC6/data/
#
# change directory to where the mapped reads are
cd $path_data/mapped_reads/ophio/
# using .bam from hisat2 mapping>samtools sort and convert 
# NO CUFFLINKS step used (nor cuffquant).
# sort and convert Hisat2 sams to sorted bams
for file in *.sam
do
samtools sort \
    -o ${file%%.*}.bam \
    ${file}
done
#
#
### Cuffdiff normalization and counting
#
# name the first sample 
# (which has 02 in the file name because first 
# sample is timepoint ZT2)
sample1=*02A*.bam
#
## Run cuffdiff for sample_TP vs all samples
for sample2 in *.bam
do
    ## run cuffdiff
        cuffdiff \
                -o ${path_data}/cuffdiff_out/ \
                -b ${path_data}/genome/ophio/ophcf2_genome.fna \
                -p 16 \
                -u ${path_data}/genome/ophio/ophcf2_genome.gff \
                ${sample1} \
                ${sample2}
#
        # rename file with gene expression
        # first make variable for pairwise name of file
        pw_name1=$(echo $sample1 | cut -d '-' -f 2 | cut -d "_" -f 1)
        #pw_name1=$(echo $sample1 | cut -d '_' -f 4,5 | cut -d '.' -f 1)
        pw_name2=$(echo $sample2 | cut -d '-' -f 2 | cut -d '_' -f 1)
        #pw_name2=$(echo $sample2 | cut -d '_' -f 4,5 | cut -d '.' -f 1)
#
   ## rename file
        cp ${path_data}/cuffdiff_out/gene_exp.diff \
		${path_data}/00_norm_gene_exp/ophio/gene_exp_${pw_name1}_vs_${pw_name2}.csv
#
## Remove all the exessive files from output folder
rm ${path_data}/cuffdiff_out/*
#
done


#
# Now run get_gene_exp_csv.py to obtain a csv with all the expression values of the samples
#
### chunk ends
