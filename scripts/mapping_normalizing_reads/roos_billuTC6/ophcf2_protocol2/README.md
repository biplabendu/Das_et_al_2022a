# Protocol 
Obtaining normalized gene expression values from short-read RNA-seq data
Used for time-course RNA-seq data of Ophiocordyceps camp-florani during manipulation

## Introduction
This documentation provides a protocol for obtaining the normalized gene expression values of time-course RNA-seq data of O. camp-florani. This involves trimming, mapping, counting, and normalization of the reads. The reads will be trimmed with BBDuk  (Bushnell B.) and mapped to the reference genome with HISAT2 (Kim et al. 2019). To obtain normalized gene expression values, we use CuffDiff (FPKM) (Trapnell et al. 2013). Each step will be carried out by a separate script. This documentation can be used to correctly run the script. However, it is recommended to read to scripts for a better understanding of the code and the analysis. The RNA-seq data is available on <insert> and a description about the RNA-seq experiment and research objectives can be found on <insert>. The genome of O. camp-florani is publicly available on <insert>. The scripts can be found on <github>.

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
mkdir  $path/TC6/data/
mkdir  $path/TC6/scripts/ # optional you can store the Ophio scripts in here
mkdir  $path/TC6/results/
```

In both the data and result directory, make a directory for the species we are analyzing, in this case O. camp-florani.
```bash
species_name = ophio
mkdir  $path/TC6/data/$species_name/
mkdir  $path/TC6/results/$species_name/
```

### Download data

```bash
mkdir /path/experiment_name/data/species_name/raw_reads
```

Also make a directory in /path/experiment_name/data/ named genome and then save the genome and annotation files in this folder.

```bash
mkdir /path/experiment_name/data/genome
```

### Trimming
This scripts excutes a adapter and quality trim of the raw reads with bbduk. The trimmed reads are stored in a trimmed_reads output directory. 
After the reads are trimmed, it will run fastqc, to check the quality of the reads. This output is stored in a fastqc_output directory.
*This scripts needs the software of bbduk and fastqc*

The used parameters are:
* "gtrim=10" will quality trim at Q10 using Phred algorithm, "qtrim=rl" will trim right and left sides
* "ref=truseq.fq.gz" included Illumina truseq adapters within BBMap package
* "hdist=1" allows one mismatch, "k=25" is k-mer size of 25, "mink=11" allow shorter k-mers of 11 at end of read
* default looks at reverse-complement and forward seuence of refseq otherwise "rcomp=<t/r>"
* "-Xmx1g" flag tells BBDuk to use 1GB of RAM
* tbo" (trim adapters based on pair overlap detection) and "tpe" (both reads to the same length) flags for normal paired-end fragment libraries