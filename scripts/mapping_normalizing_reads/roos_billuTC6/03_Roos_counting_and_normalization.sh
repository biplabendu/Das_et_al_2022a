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

# Set directory
cd $path/mapped_reads/ophio/

#
#nested loop for pairwise comparison of TC2 vs TC4, TC2 vs TC6 ect., TC4 vs TC6, TC4 vs TC8, ect.
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

#
# TO DO: Check HTSeq-based TMM normalization instead of Cuffdiff-based FPKM/RPKM normalization
#
### chunk ends