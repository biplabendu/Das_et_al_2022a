# Housekeeping ---------------------------------------------------------------
#
set.seed(420)
rm(list = ls())
#
## Load packages ----------
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
#
# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
#
## set parameters and thresholds --------
#
# gamma-pvalue threshold for inferring rhythmicity
gamma.pval = 0.05
#
# path to repo
path = "/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a"

## Path to save files ------
# for supplementary files
# supp.path="~/University\ of\ Central\ Florida/Charissa\ De\ Bekker\ -\ Ant-Fungus-Clock-Interactions/04_manuscript/03_supplementary_files/"

# 01. Databases -----------------------------------------------------------

# 1. TC6_fungal_ejtk.db
#
# Desc: This database will contain all ejtk-output for fungal expression data collected for TC6
#
### Contents (last updated 20-Sep-21)
# a. ophio_cflo_zscores_24h
# b. ophio_cflo_zscores_12h
# c. ophio_cflo_zscores_08h
# d. ophio_kim_LD_zscores_24h
# e. ophio_kim_DD_zscores_24h
#
# Load the data %
my.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC6_fungal_ejtk.db")
# which tables are in the database
#
src_dbi(my.db)

# Load library glue for naming purposes
library(glue)
#
# Choose species
species <- 'beau'
## End.

## z-score >> eJTK ------------------------
# 
# Give the periods (8h, 12h, 24h)
periods <- c('08',12,24) %>%
  as.character

for (i in periods) {
  # name the file with the species and period
  name <- glue('{species}_zscores_{i}h')
  
  # CSV name
  csv.name <- glue("./results/ejtk_output/{species}/zscore/{species}_zscores_noNAs_cos{i}_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt")
  #
  # Note: ophio.cflo.24.zscore is exactly the same as ophio.cflo.24 meaning that using
  #       FPKM values or z-scores to run eJTK does not change the result

  # Period
  zscore <- read.csv(csv.name,
                     sep = "\t", header = T, stringsAsFactors = F)
  zscore %>%
    head()
  # Save file to database
  dbWriteTable(my.db, name, zscore)
}

# Rhytmic or not ----------------------------------------------------------
# tbl(my.db, “tbl_name”) %>% collect() %>% view() for viewing table
#
### Write table with all the gene ID's in it and say if it is rhytmic or not
period <- '24'

# path to file with all the gene IDs
csv.name <- paste0('./data/input/', species, '/', species, '_gene_names_robin_ncbi.csv')

# name of the new table in the database
name <- glue('{species}_rhythmic_genes_{period}h')

# read in the gene IDs
gene.IDs <- read.csv(csv.name, sep= ',', header = TRUE, stringsAsFactors = FALSE)

# Make the table
# For rhythmic genes
#
# each species needs per period a table, for instance named beu_rhytmic_genes_24h ect.
# Tables should have 3 cols; RobinID NcbiID, rhytmic == yes/no

rhythmic <- 
  my.db %>% 
  tbl(glue("{species}_zscores_{period}h")) %>% 
  select(gene_ID_NCBI = ID, 
         GammaP) %>% 
  collect() %>% 
  
  # mutate(my_col_name = GammaP*2) %>% you can do cool stuff with mutate
  
  mutate(rhythmic = ifelse(GammaP < 0.05, "yes", "no")) %>%   
  
  select(-GammaP) %>%

  left_join(gene.IDs, by = c('gene_ID_NCBI' = "attributes_ncbi")) %>% 
  
  select(-c(start,end))

# group_by(rhythmic_24h) %>% 
# 
# summarise(n_genes = n())
#### 
#

# Make the table in the database
dbWriteTable(my.db, name, rhythmic)

##### Rhytmicity for ophio_kim ----------------------

gene.ID <- read_csv(glue('{path}/data/input/ophio_kim/ophio_kim_annots_robin_ncbi.csv')) %>% 
  select(gene_ID_robin, gene_ID_ncbi)

rhythmic <- 
  my.db %>% 
  tbl('ophio_kim_LD_zscores_24h') %>% 
  select(gene_ID_ncbi = ID, 
         GammaP) %>% 
  collect() %>% 
  mutate(rhythmic = ifelse(GammaP < 0.05, "yes", "no")) %>% 
  left_join(gene.ID, by = 'gene_ID_ncbi')

dbWriteTable(my.db, 'ophio_kim_LD_rhythmic_genes_24h', rhythmic)

##### Done.


