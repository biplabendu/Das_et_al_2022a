#!/bin/bash
#SBATCH -J mapping
#SBATCH -t 4:0:0
#SBATCH --mem=5G


#
# fresh ophio index, with exons
gffread "ophcf2_genome.gff" -T -o ophcf2.gtf
hisat2_extract_splice_sites.py "ophcf2.gtf" > ophcf2_splicesites.txt
hisat2_extract_exons.py "ophcf2.gtf" > ophcf2_exons.txt
#
# echo "BUILDING THE OPHIO INDEX WITH SPLICE AND EXONS SITE"
# index the exon site version of the ophcf2 genome
hisat2-build -f --ss ophcf2_splicesites.txt --exon ophcf2_exons.txt ophcf2_genome.fna ophcf2_exons_index
#
# map (trimmed) rnaseq sample to (indexed) genome
#
### TO DO: MAKE A LOOP TO INDIVIDUALLY MAP ALL RNASEQ SAMPLES TO THE GENOME
#
mkdir $path/mapped_reads/ophio/
#
cd $path/trimmed_reads/ophio/
for file in *fastq.gz
do
hisat2 \
        -p 8 \
        --new-summary \
        -q \
        --dta-cufflinks \
        -x $path/genome/ophio/ophcf2_exons_index \
        -U $path/trimmed_reads/ophio/${file} \
        -S $path/mapped_reads/ophio/mapped_${file%%.*}.sam
done
#
### chunk ends
