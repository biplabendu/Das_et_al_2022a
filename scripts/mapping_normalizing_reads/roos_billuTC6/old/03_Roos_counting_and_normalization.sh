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

## Sort and convert Hisat2 obtained sams to sorted bams
for file in $path/mapped_reads/ophio/*.sam
do
samtools sort -o ${file%%.*}.bam ${file}
done

### Cuffdiff normalization and counting
#
# Make directoty for output cuffdiff
mkdir $path/cuffdiff_out
mkdir $path/cuffdiff_out/ophio
#
# name the first sample (which has 2 in the file name because first sample is timepoint 2)
sample1=mapped_trimmed_AQ_*02*.bam
#
## Run cuffdiff for sample_TP vs all samples
for sample2 in *.bam
do
        cuffdiff \
                -o ${path}/cuffdiff_out/ophio/ \
                -b ${path}/genome/ophio/ophcf2_genome.fna \
                -p 16 \
                -u ${path}/genome/ophio/ophcf2_genome.gff \
                ${sample1} \
                ${sample2}

        # rename file with gene expression
        # first make variable for pairwise name of file
        pw_name1=$(echo $sample1 | cut -d '-' -f 2 | cut -d "_" -f 1)
        pw_name2=$(echo "$sample2" | cut -d '-' -f 2 | cut -d '_' -f 1)
        # rename file
        mv ${path}/cuffdiff_out/ophio/gene_exp.diff ${path}/cuffdiff_out/ophio/gene_exp_${pw_name1}vs${pw_name2}.csv
done

## Remove all the exessive files from output folder
cd ${path}/cuffdiff_out/ophio/
rm !(*.csv)   
#
# Now run get_gene_exp_csv.py to obtain a csv with all the expression values of the samples
#
### chunk ends