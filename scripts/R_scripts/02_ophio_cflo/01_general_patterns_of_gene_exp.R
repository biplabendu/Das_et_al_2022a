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
## load functions ---------
source("./functions/enrichment_analysis.R")
#
## set parameters and thresholds --------
#
# gamma-pvalue threshold for inferring rhythmicity
gamma.pval = 0.05
#
## Path to save files
# for supplementary files
supp.path="~/University\ of\ Central\ Florida/Charissa\ De\ Bekker\ -\ Ant-Fungus-Clock-Interactions/04_manuscript/03_supplementary_files/"

# 00. Databases -----------------------------------------------------------
#
# 1. TC6_fungal_ejtk.db
# Desc: This database will contain all ejtk-output for fungal expression data collected for TC6
ejtk.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC6_fungal_ejtk.db")
# which tables are in the database
src_dbi(ejtk.db)
#
# 2. TC6_fungal_data.db
data.db <- dbConnect(RSQLite::SQLite(),
                     "./data/databases/TC6_fungal_data.db")
src_dbi(data.db)
#
##


# 01. General patterns of gene expression ---------------------------------

# A1: genes that have NO expression (FPKM == 0 at all time points)
not.expressed <-
  tbl(data.db, "ophio_cflo_fpkm") %>% 
  collect() %>% #nrow()
  filter_at(vars(starts_with("Z")), all_vars(. == 0)) %>% 
  pull(gene_name)
# A2: run enrichment
not.expressed %>% 
  go_enrichment(., 
                org = "ophio_cflo", 
                bg = "all") # enrichment against all ophio_cflo genes in the genome

# B: genes that are expressed (FPKM > 1 for at least one time point)
expressed <- 
  tbl(data.db, "ophio_cflo_expressed_genes") %>% 
  filter(expressed=="yes") %>% 
  collect() %>% 
  pull(gene_name)


# 02. Diurnal rhythms in gene expression ----------------------------------

# A1. number of 24h-rhythmic genes at GammaP < 0.05
rhy.24 <- tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% filter(GammaP < gamma.pval) %>% collect()

# A2. save the ejtk-results (ALL GENES) as a supplementary file
tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% 
  collect() %>% 
  write.csv(.,
            file = paste0(supp.path,"S1_ejtk_24h_ophio_cflo_all_genes.csv"),
            row.names = F)

# A3. save the ejtk-results (rhy24 genes only) as a supplementary file
tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% 
  filter(GammaP < gamma.pval) %>% 
  collect() %>% 
  write.csv(.,
            file = paste0(supp.path,"S2_ejtk_24h_ophio_cflo_rhy_genes.csv"),
            row.names = F)
  