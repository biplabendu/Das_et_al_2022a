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
# Load the data
my.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC6_fungal_ejtk.db")
# which tables are in the database
#
src_dbi(my.db)

# Load library glue for naming purposes
library(glue)
#
# Choose species
species <- 'ophio_cflo'
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

