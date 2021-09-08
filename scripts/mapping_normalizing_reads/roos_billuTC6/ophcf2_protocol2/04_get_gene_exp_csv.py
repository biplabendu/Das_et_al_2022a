#! /usr/bin/env python3
#SBATCH -J get_csv
#SBATCH -t 4:0:0
#SBATCH --mem=5G

### Get one csv with the normalized expression data 
# This Python script to make a csv of the output data from cuffdiff. It extracts the gene_id and expressionvalues for each analysis and writes it to 1 csv file
# This scripts needs os.system and pandas

# 00: import system commands and pandas (for DF management)
import os
import pandas as pd
#
# set path were data is
path = '/home/uu_bio_fg/rbrouns/data/ant_fungus/TC6/data/ophio/cuffdiff_out/'

# 01: Create a list with all the files in the directory where the gene_ex data is
ls_gene_ex_files = (os.listdir(path))
# sort the list
ls_gene_ex_files.sort()
# some how the sample for file 20 does not contain data so that one is removed from te list
ls_gene_ex_files.remove('gene_exp_02Avs20A.csv')

### 02: Create DF frame with rownames and without exp_values
#
# get first filename of list (thus file of sample 1)
sample1 = (ls_gene_ex_files[0])
# create df with the data from sample 1
file0 = pd.read_csv(f'{path}{sample1}', sep='\t')
# extract the columns test_id, gene_id, gene, locus into new DF
TC6_gene_ex = file0[['gene_id', 'gene', 'locus']]

### 03: Append expression values as new column matchin to rownames
#
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

# 04: Write DF to csv file
TC6_gene_ex.to_csv(f'{path}TC6_gene_exp.csv')
#
### End
#
### To run in command line ./script_name.py 