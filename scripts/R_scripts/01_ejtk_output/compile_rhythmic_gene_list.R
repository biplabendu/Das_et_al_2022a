
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
supp.path="~/University\ of\ Central\ Florida/Charissa\ De\ Bekker\ -\ Ant-Fungus-Clock-Interactions/04_manuscript/03_supplementary_files/"


# 01. Databases -----------------------------------------------------------

# 1. TC6_fungal_ejtk.db
#
# Desc: This database will contain all ejtk-output for fungal expression data collected for TC6
#
### Contents (last updated 02-Sep-21)
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
src_dbi(my.db)
#
## End.

# 02. Ophiocordyceps camponoti-floridani --------------------------------------
#
## FPKM >> eJTK ------------------------
# #
# #
# ## Period = 24h
# ophio.cflo.24 <- read.csv("./results/ejtk_output/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
#            sep = "\t", header = T, stringsAsFactors = F)
# ophio.cflo.24 <- ophio.cflo.24 %>%  
#   filter(GammaP < gamma.pval) %>% 
#   pull(ID)
# #
# #
# ## Period = 12h
# ophio.cflo.12 <- read.csv("./results/ejtk_output/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples_cos12_ph0020by4_as0420by4_jtkout_GammaP.txt",
#                           sep = "\t", header = T, stringsAsFactors = F)
# ophio.cflo.12 <- ophio.cflo.12 %>% 
#   filter(GammaP < gamma.pval) %>% 
#   pull(ID)
# #
# #
# ## Period = 8h
# ophio.cflo.08 <- read.csv("./results/ejtk_output/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples_cos8_ph0020by4_as0420by4_jtkout_GammaP.txt",
#                           sep = "\t", header = T, stringsAsFactors = F)
# ophio.cflo.08 <- ophio.cflo.08 %>% 
#   filter(GammaP < gamma.pval) %>% 
#   pull(ID)
# #
# #
# ## Save the rhythmic gene-lists as a data file
# save(ophio.cflo.08,
#      ophio.cflo.12,
#      ophio.cflo.24,
#      file = "./data/rdata/rhythmic_genes/ophio_cflo_rhy_genes_gammap_05.RData")
# #
# ## To read the RData file, use:
# ## load("/path/to/file.RData)
# ## End.
# #
# #
## z-score >> eJTK ------------------------
#
# Note: ophio.cflo.24.zscore is exactly the same as ophio.cflo.24 meaning that using
#       FPKM values or z-scores to run eJTK does not change the result
#
## Period = 24h
ophio.cflo.24.zscore <- read.csv("./results/ejtk_output/ophio_cflo/from_zscore/ophio_cflo_zscores_noNAs_cos24_ph0020by4_as0420by4_zscore_jtkout_GammaP.txt",
                          sep = "\t", header = T, stringsAsFactors = F)
ophio.cflo.24.zscore %>% 
  head()
# Save file to database
dbWriteTable(my.db, "ophio_cflo_zscores_24h", ophio.cflo.24.zscore)
#
#
## Period = 12h
ophio.cflo.12.zscore <- read.csv("./results/ejtk_output/ophio_cflo/from_zscore/ophio_cflo_zscores_noNAs_cos12_ph0020by4_as0420by4_zscore_jtkout_GammaP.txt",
                          sep = "\t", header = T, stringsAsFactors = F)
ophio.cflo.12.zscore %>% 
  head()
# Save file to database
dbWriteTable(my.db, "ophio_cflo_zscores_12h", ophio.cflo.12.zscore)
#
#
## Period = 8h
ophio.cflo.08.zscore <- read.csv("./results/ejtk_output/ophio_cflo/from_zscore/ophio_cflo_zscores_noNAs_cos08_ph0020by4_as0420by4_zscore_jtkout_GammaP.txt",
                          sep = "\t", header = T, stringsAsFactors = F)
ophio.cflo.08.zscore %>% 
  head()
# Save file to database
dbWriteTable(my.db, "ophio_cflo_zscores_08h", ophio.cflo.08.zscore)
#
#
## End.


# 03. Ophiocordyceps kimflemingae ---------------------------------------------
#
## FPKM >> eJTK ------------------------
# #
# ## Light-Dark cycle
#   ## Period = 24h
# ophio.kim.ld.24 <- read.csv("./results/ejtk_output/ophio_kim/ld/normalized_gene_exp_ophio_kim_ld_samples_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
#                             sep = "\t", header = T, stringsAsFactors = F)
# ophio.kim.ld.24 <- ophio.kim.ld.24 %>% 
#   filter(GammaP < gamma.pval) %>% 
#   pull(ID)
# #
# #
# ## Dark-Dark
#   ## Period = 24h
# ophio.kim.dd.24 <- read.csv("./results/ejtk_output/ophio_kim/dd/normalized_gene_exp_ophio_kim_dd_samples_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
#                             sep = "\t", header = T, stringsAsFactors = F)
# ophio.kim.dd.24 <- ophio.kim.dd.24 %>% 
#   filter(GammaP < gamma.pval) %>% 
#   pull(ID)
# #
# #
# ## Save the rhythmic gene-lists as a data file
# save(ophio.kim.dd.24,
#      ophio.kim.ld.24,
#      file = "./data/rdata/rhythmic_genes/ophio_kim_rhy_genes_gammap_05.RData")
# #
# ##
## End.
#
#
## z-score >> eJTK ------------------------
#
## Light-Dark cycle
  ## Period = 24h
ophio.kim.ld.24.zscore <- read.csv("./results/ejtk_output/ophio_kim/ld/from_zscore/ophio_kim_LD_zscores_noNAs_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
ophio.kim.ld.24.zscore %>% 
  mutate(experiment="light-dark") %>% 
# save file to database
  dbWriteTable(my.db, "ophio_kim_LD_zscores_24h", .)
#
#
## Dark-Dark
  ## Period = 24h
ophio.kim.dd.24.zscore <- read.csv("./results/ejtk_output/ophio_kim/dd/from_zscore/ophio_kim_DD_zscores_noNAs_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
ophio.kim.dd.24.zscore %>% 
  mutate(experiment="dark-dark") %>%
# save file to database
  dbWriteTable(my.db, "ophio_kim_DD_zscores_24h", .)
#
##
## End.
