
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
## set parameters and thresholds
#


# 00. Databases -----------------------------------------------------------

# 1. TC6_fungal_data.db
# Desc: This database will contain all fungal expression data collected for TC6
### Contents -
# a. Ophio_cflo raw fpkm
# b. Ophio_cflo raw fpkm + expressed (yes/no)
# c. 
#
# Load the data
my.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC6_fungal_data.db")
# which tables are in the database
src_dbi(my.db)

# 01. Ophio_cflo ----------------------------------------------------------

o.cflo <- read.csv("./results/normalized_gene_exp/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples.csv",
                   header = T, stringsAsFactors = F, na.strings = c(NA, "", " "))
# # Save file to database
# dbWriteTable(my.db, "ophio_cflo_fpkm", o.cflo)


# Expressed genes ---------------------------------------------------------
# list of all ophio_cflo "Expressed" genes (>1 FPKM during the 24h period)
expressed <-
  o.cflo %>% 
  na.omit() %>%  
  filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>% 
  pull(gene_name) %>% 
  unique()
# append this information to the o.cflo data
expressed.o.cflo <-
  o.cflo %>% 
  select(gene_name) %>% 
  mutate(expressed = ifelse(gene_name %in% expressed,"yes","no"))
# # Save file to databse
# dbWriteTable(my.db, "ophio_cflo_expressed_genes", expressed.o.cflo)

# log2-expression --------------------------------------------------------
# log2-transform the data
gene.names <- o.cflo %>% pull(gene_name)
log2.o.cflo <- log2(o.cflo[-1] + 1)
log2.o.cflo$gene_name <- gene.names
log2.o.cflo <- log2.o.cflo %>% select(gene_name, everything())
# check the log2-transformed data
log2.o.cflo %>% str()
# # Save file to database
# dbWriteTable(my.db, "ophio_cflo_log2fpkm", log2.o.cflo)

# zscore-log2-expression --------------------------------------------------
# z-score the data
gene.names <- log2.o.cflo %>% pull(gene_name)
sample.names <- names(log2.o.cflo[-1])
zscores.o.cflo <- log2.o.cflo %>%
  # create a gene x exp matrix
  select(-1) %>% 
  as.matrix() %>%
  # use the scale function for each row to calculate z-scores
    # scale calculates (x-mean(X))/sd(X)
    # 1 indicates row-wise
  apply(., 1, scale) %>% 
  # the output needs to be transposed
  t() %>% 
  # make it a dataframe 
  as.data.frame()
# add the column names
names(zscores.o.cflo) <- sample.names
zscores.o.cflo$gene_name <- gene.names
zscores.o.cflo <- zscores.o.cflo %>% select(gene_name, everything())
# check the z-score transformed dataset
zscores.o.cflo %>% str()
# # Save file to database
# dbWriteTable(my.db, "ophio_cflo_zscores", zscores.o.cflo)
