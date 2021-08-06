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
path=/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/mapped_reads/ophio

# using .bam from hisat2 mapping>samtools sort and convert 
# NO CUFFLINKS step used (nor cuffquant).
# sort and convert Hisat2 sams to sorted bams
for file in $path/*.sam
do
samtools sort \
    -o ${file%%.*}.bam \
    ${file}
done

#
# use the following code to obtain normalized FPKM values
# TO DO: Make a loop to obtain normalized gene counts for each sample
#
cuffdiff -o Foragers -b "/home/billu/cflo_genome/camp_genome.fna" -p 16 -u "/home/billu/cflo_genome/GCF_003227725.1_Cflo_v7.5_genomic.gff" \2F_time5_1.bam \4F_time5_1.bam
#
# TO DO: Check HTSeq-based TMM normalization instead of Cuffdiff-based FPKM/RPKM normalization
#
### chunk ends