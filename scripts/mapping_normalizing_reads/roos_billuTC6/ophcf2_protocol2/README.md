# Protocol 
Obtaining normalized gene expression values from short-read RNA-seq data
Used for time-course RNA-seq data of Ophiocordyceps camp-florani during manipulation

## Introduction
This documentation provides a protocol for obtaining the normalized gene expression values of time-course RNA-seq data of O. camp-florani. This involves trimming, mapping, counting, and normalization of the reads. The reads will be trimmed with BBDuk  (Bushnell B.) and mapped to the reference genome with HISAT2 (Kim et al. 2019). To obtain normalized gene expression values, we use CuffDiff (FPKM) (Trapnell et al. 2013). This documentation provides a step by step explanation of all the code and can used to run the analysis. Each step can also be carried out by separate scripts, which are available on the git. The RNA-seq data is available on <insert> and a description about the RNA-seq experiment and research objectives can be found on <insert>. The genome of O. camp-florani is publicly available on <insert>. The scripts can be found on <github>.

*References:*
*BBMap – Bushnell B. – sourceforge.net/projects/bbmap/*
*Kim, Daehwan, Joseph M Paggi, Chanhee Park, Christopher Bennett, and Steven L Salzberg. 2019. “Graph-Based Genome Alignment and Genotyping with HISAT2 and HISAT-Genotype.” Nature Biotechnology 37 (8): 907–15.*
*Trapnell, Cole, David G Hendrickson, Martin Sauvageau, Loyal Goff, John L Rinn, and Lior Pachter. 2013. “Differential Analysis of Gene Regulation at Transcript Resolution with RNA-Seq.” Nature Biotechnology 31 (1): 46–53.*


## Requirements
This analysis requires the following software/modules:
* BBMap
* HISAT2
* FastQC
* Samtools
* CuffLinks
* python 3.7 with modules os and pandas

This analysis requires the following scripts:
* 00_data_strusture.sh
* 00_data_download.sh
* 01_trim_reads.sh
* 02_mapping.sh
* 03_counting_and_normalization.sh
* 04_get_gene_exp_csv.py


## Usage

### Data structure
The 00_data_structure.sh script creates the directories to organize the data. It is important to for executing the other scripts.

Here we defined the absolute path to where the directory needs to be created. 

```bash
path = /home/uu_bio_fg/rbrouns/data/ant_fungus # don’t add a / at the end of the path

```

Then create a directory for the experiment,in this case Time-course experiment 6 (TC6)
```bash
experiment_name = TC6
mkdir $path/$experiment_name/
```

In the experiment directory, create a directory for storing the data, scripts and the results
```bash
mkdir  $path/$experiment_name/data/
mkdir  $path/$experiment_name/scripts/ # optional you can store the Ophio scripts in here
mkdir  $path/$experiment_name/results/
```

In both the data and result directory, make a directory for the species we are analyzing, in this case O. camp-florani.
```bash
species_name = ophio
mkdir  $path/TC6/data/$species_name/
mkdir  $path/TC6/results/$species_name/
```

### Download data
The next step is to download the RNA seq data, the genome file and annotation file. The RNA seq is available on ... The genome files are publicity available on ...

First we creata a folder to store the RNA seq reads.

```bash
mkdir $path/$experiment_name/data/$species_name/raw_reads
```

Then, this commands downloads the RNA seq data into the right folder
```bash
<ad code here>
```

Additionally, make a folder in /path/experiment_name/data/ named genome to save the genome and annotation files in.

```bash
mkdir $path/$experiment_name/data/$species_name/genome
cd $path/$experiment_name/data/$species_name/genome/
```

The next command downloads the genome file of Ophio camp-f from the NCBI site 
```bash
# download the genome file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.fna.gz
# unzip the genome file
gzip -d GCA_012980515.1_Ophcf2_genomic.fna.gz
# rename the genome file 
mv GCA_012980515.1_Ophcf2_genomic.fna ophcf2_genome.fna
```

Next, download the annotation file, unzip and rename it.
```bash
# download
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/980/515/GCA_012980515.1_Ophcf2/GCA_012980515.1_Ophcf2_genomic.gff.gz
# unzip 
gzip -d GCA_012980515.1_Ophcf2_genomic.gff.gz
# rename 
mv GCA_012980515.1_Ophcf2_genomic.gff ophcf2_genome.gff
```

### Trimming
After the data is downloaded, we can start with trimming the data. This is done to remove left over adapter sequences and improve the quality of the reads. This scripts excutes a adapter and quality trim of the raw reads with bbduk. The trimmed reads are stored in a trimmed_reads output directory. 
After the reads are trimmed, it will run fastqc, to check the quality of the reads. This output is stored in a fastqc_output directory.
*This scripts needs the software of bbduk and fastqc*

First, make output folder for the trimmed reads and set as working directory
```bash
mkdir $path/$experiment_name/data/ophio/trimmed_reads
cd /$path/$experiment_name/data/ophio/demultiplexed_reads/
```

Then do the trimming of all reads in one go.
For this bit of code the file with the used adapeters is needed at the parameter ref=<file.fa>. This fasta file contains the sequences of the used adapters during sequencing. (The path used in this script is al follows, /home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/bbmap_adapters/adapters.fa)

First we do the adapter trimming, then the quality trimming
```bash
for file in *fastq.gz
do
# adapter trimming
bbduk.sh \
	-Xmx1g \
	in1=${file} \
	out1=./../ophio/trimmed_reads/trimmed_A_${file} \
	ref=/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/bbmap_adapters/adapters.fa \
	ktrim=r \
	k=23 \
	mink=11 \
	hdist=1 \
	tpe \
	tbo	
# quality trimming
bbduk.sh \
	-Xmx1g \
	in=./../ophio/trimmed_reads/trimmed_A_${file} \
	out=./../ophio/trimmed_reads/trimmed_AQ_${file} \
	qtrim=rl \
	trimq=10
done
```
The used parameters are:
* "gtrim=10" will quality trim at Q10 using Phred algorithm, "qtrim=rl" will trim right and left sides
* "ref=truseq.fq.gz" included Illumina truseq adapters within BBMap package
* "hdist=1" allows one mismatch, "k=25" is k-mer size of 25, "mink=11" allow shorter k-mers of 11 at end of read
* default looks at reverse-complement and forward seuence of refseq otherwise "rcomp=<t/r>"
* "-Xmx1g" flag tells BBDuk to use 1GB of RAM
* tbo" (trim adapters based on pair overlap detection) and "tpe" (both reads to the same length) flags for normal paired-end fragment libraries

Next remove all the trimmed_A_* files, which are the files that are only trimmed on adapters. What remains are the adapter and quality trimmed reads.
```bash
rm ./../trimmed_reads/ophio/trimmed_A_*
```
### Mapping 
After the reads are trimmed, the next step is mapping the reads to the refenrence genome. This script maps the reads to the reference genome, with HISAT2 (ref). The refence genomes are publicity available on ...
*For this script the software of HISAT2 and gffread is needed.*

We will first start with indexing the O. camp-florani genome, with exon and and splice sites. 

First, set genome folder as working directory to store the index files we're creating
```bash
cd $path/$experiment_name/data/$species_name/genome
```
Since we work with RNA-seq data (and thus only transcibed DNA, which are the introns), we want to index the exons and slice sites for the mapping.
```bash
gffread "ophcf2_genome.gff" -T -o ophcf2.gtf 
hisat2_extract_splice_sites.py "ophcf2.gtf" > ophcf2_splicesites.txt
hisat2_extract_exons.py "ophcf2.gtf" > ophcf2_exons.txt
```
Then we index the exon site version of the ophcf2 genome with hisat2-build
```bash
hisat2-build \ 
    -f \ # The reference input files are FASTA files
    --ss ophcf2_splicesites.txt \ # Input the splice sites
    --exon ophcf2_exons.txt \ # Input the exons
    ophcf2_genome.fna \ # Input the reference genome
    ophcf2_exons_index # Basename of output file to write indexed genome with exons and splice sites in
```

The next step is mapping the trimmed RNAseq sample reads onto the indexed genome
```bash
# Make a folder to store the mapped reads
mkdir $path/$experiment_name/data/ophio/mapped_reads
# And set the folder with trimmed reads as working directory
cd $path/$experiment_name/data/ophio/trimmed_reads/
```
The next code maps the reads with hisat2 by looping over each file in the folder
```bash
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
```
### Counting and Normalization
The next step is counting the read that are mapped to the genes and normalize them. This will result in the final gene expression values. This script counts the mapped reads and normalizes them to get the gene expression values in FRKPM with CuffDiff (ref).
*This script requires Samtools and cuffdiff (from CuffLinks)*

First, sort and convert all Hisat2 sams files to sorted bams with samtools sort
```bash
# move into folder where the mapped reads are
cd $path/data/ophio/mapped_reads/

# convert sam files to bam files
for file in *.sam
do
samtools sort \
    -o ${file%%.*}.bam \ # Output same basename with .bam extension
    ${file} # Input .sam files
done
```
Next step is the Cuffdiff normalization. Cuffdiff compares one sample to another to obtain normalized gene expression values. Therefore we define our first sample and compare all the other samples to this sample. Per sample vs sample analysis cuffdiff outputs a file with the gene expressoin values. The orinalal file name of this file is gene_exp.diff, but we will rename this accoringly to which samples are compared.
```bash
# Create output folder
mkdir $path/data/ophio/cuffdiff_out

# Set working directory to mapped reads
cd ${path}/data/ophio/mapped_reads/

# Name the first sample (which has 2 in the file name because first sample is timepoint 2)
sample1=mapped_trimmed_AQ_*02*.bam

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
```
Then we will remove all the other outputs files cuffdiff creates, since were only interested in the gene expression values
```bash
# remove all the other files 
cd ${path}/data/ophio/cuffdiff_out/
rm !(*.csv) 
```

Now all the gene eprression value files will be in the cuffdiff_out folder.

### Get one csv with the normalized expression data
The next script is in Python 3.7, and needs to be ran seprately from the rest. This script makes one csv file out of all the output data from cuffdiff. It extracts the gene_id and expression values for each sample vs sample analysis. The resulting data file can be easily used for different anaysis in python or R.
*This scripts run in Python 3.7 and need the os.system and pandas modules*

00: Import system commands and pandas (for DataFrame management). The os.system modules allows us to use bash command in python. The dataframes are a part of the pandas modules, which makes data easy to handle.
```python
import os
import pandas as pd
```
Then we need to specify the path were data is (the cuddfiff_out folder)
```python
path = '/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/ophio/cuffdiff_out/'
```

01: Create a list with all the files in the directory where the gene_ex data is, so that we can retreive the data
```python
# create the list
ls_gene_ex_files = (os.listdir(path))
# sort the list
ls_gene_ex_files.sort()
# some how the sample for file 20 does not contain data so that one is removed from te list
ls_gene_ex_files.remove('gene_exp_02Avs20A.csv')
```
02: Create a empty  DataFrame (meaning no expression value data is added yet) with the rows being diffrent genes and the columns containing the gene_id, gene name and the locus
```python
# get first filename of list (thus file of sample 1)
sample1 = (ls_gene_ex_files[0])
# create df with the data from sample 1
file0 = pd.read_csv(f'{path}{sample1}', sep='\t')
# extract the columns test_id, gene_id, gene, locus into new DF
TC6_gene_ex = file0[['gene_id', 'gene', 'locus']]
```
03: Add expression values for each sample as new column in the DataFrame
```python
for sample in ls_gene_ex_files:
   # create df with the data from sample
   file = pd.read_csv(f'{path}{sample}', sep='\t')
   # extract the column of expression value_2 and 
   value_2 = file[['gene_id','value_2']]
   # rename header of expression value to sample name
   s_header = sample.split('.')[0]
   ss_header = 'sample_' + s_header.split('vs')[1] 
   value_2.columns = ['gene_id',ss_header]
   # append TC6 DF with expression values DF
   TC6_gene_ex = pd.concat([TC6_gene_ex, value_2.iloc[:,1]], axis=1)
```
04: Then write teh DataFrame to a csv file
```python
TC6_gene_ex.to_csv(f'{path}TC6_gene_exp.csv')
```

### End
Now you will have a csv file containing al the expression values of the Ophio camp-f RNA seq experiment. Each row is a different gene and each colums is a different sample, thus a different time point. This data is used to find rhytmicity in genes and diffrentially expressed genes.

