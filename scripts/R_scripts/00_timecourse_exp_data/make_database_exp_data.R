
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
#
# Desc: This database will contain all fungal expression data collected for TC6
#
### Contents (last updated: 02-Sep-2021)
# a. ophio_cflo_expressed_genes
# b. ophio_cflo_fpkm
# c. ophio_cflo_log2fpkm
# d. ophio_cflo_zscores
# e. ophio_kim_DD_expressed_genes
# f. ophio_kim_DD_fpkm
# g. ophio_kim_DD_log2fpkm
# h. ophio_kim_DD_zscores
# i. ophio_kim_expressed_genes
# j. ophio_kim_fpkm
# k. ophio_kim_log2fpkm
# l. ophio_kim_zscores
#
# Load the data
my.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC6_fungal_data.db")
# which tables are in the database
src_dbi(my.db)

# # 01. Ophio_cflo ----------------------------------------------------------
# 
# o.cflo <- read.csv("./results/normalized_gene_exp/raw_fpkm/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples.csv",
#                    header = T, stringsAsFactors = F, na.strings = c(NA, "", " "))
# # # Save file to database
# # dbWriteTable(my.db, "ophio_cflo_fpkm", o.cflo)
# 
# 
# # Expressed genes ---------------------------------------------------------
# # list of all ophio_cflo "Expressed" genes (>1 FPKM during the 24h period)
# expressed <-
#   o.cflo %>% 
#   na.omit() %>%  
#   filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>% 
#   pull(gene_name) %>% 
#   unique()
# # append this information to the o.cflo data
# expressed.o.cflo <-
#   o.cflo %>% 
#   select(gene_name) %>% 
#   mutate(expressed = ifelse(gene_name %in% expressed,"yes","no"))
# # # Save file to databse
# # dbWriteTable(my.db, "ophio_cflo_expressed_genes", expressed.o.cflo)
# 
# # log2-expression --------------------------------------------------------
# # log2-transform the data
# gene.names <- o.cflo %>% pull(gene_name)
# log2.o.cflo <- log2(o.cflo[-1] + 1)
# log2.o.cflo$gene_name <- gene.names
# log2.o.cflo <- log2.o.cflo %>% select(gene_name, everything())
# # check the log2-transformed data
# log2.o.cflo %>% str()
# # # Save file to database
# # dbWriteTable(my.db, "ophio_cflo_log2fpkm", log2.o.cflo)
# 
# # zscore-log2-expression --------------------------------------------------
# # z-score the data
# gene.names <- log2.o.cflo %>% pull(gene_name)
# sample.names <- names(log2.o.cflo[-1])
# zscores.o.cflo <- log2.o.cflo %>%
#   # create a gene x exp matrix
#   select(-1) %>% 
#   as.matrix() %>%
#   # use the scale function for each row to calculate z-scores
#     # scale calculates (x-mean(X))/sd(X)
#     # 1 indicates row-wise
#   apply(., 1, scale) %>% 
#   # the output needs to be transposed
#   t() %>% 
#   # make it a dataframe 
#   as.data.frame()
# # add the column names
# names(zscores.o.cflo) <- sample.names
# zscores.o.cflo$gene_name <- gene.names
# zscores.o.cflo <- zscores.o.cflo %>% select(gene_name, everything())
# # check the z-score transformed dataset
# zscores.o.cflo %>% str()
# # # Save file to database
# # dbWriteTable(my.db, "ophio_cflo_zscores", zscores.o.cflo)

# # Save a csv with the zscores
# ocflo.zscore <- tbl(my.db, "ophio_cflo_zscores") %>% collect()
# 
# ocflo.zscore.noNAs <- 
#   ocflo.zscore %>% 
#   na.omit()
# 
# write.csv(ocflo.zscore.noNAs,
#           file = "./results/normalized_gene_exp/zscore/ophio_cflo/ophio_cflo_zscores_noNas.csv",
#           row.names = F)


# 02. Ophio_kim LD -----------------------------------------------------------

# o.kim <- read.csv("./results/normalized_gene_exp/raw_fpkm/ophio_kim/ld/normalized_gene_exp_ophio_kim_ld_samples.csv",
#                    header = T, stringsAsFactors = F, na.strings = c(NA, "", " "))
# # # Save file to database
# # dbWriteTable(my.db, "ophio_kim_fpkm", o.kim)
# 
# 
# # Expressed genes ---------------------------------------------------------
# # list of all ophio_cflo "Expressed" genes (>1 FPKM during the 24h period)
# expressed <-
#   o.kim %>%
#   na.omit() %>%
#   filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>%
#   pull(gene_name) %>%
#   unique()
# # append this information to the o.kim data
# expressed.o.kim <-
#   o.kim %>%
#   select(gene_name) %>%
#   mutate(expressed = ifelse(gene_name %in% expressed,"yes","no"))
# # # Save file to databse
# # dbWriteTable(my.db, "ophio_kim_expressed_genes", expressed.o.kim)
# 
# # log2-expression --------------------------------------------------------
# # log2-transform the data
# gene.names <- o.kim %>% pull(gene_name)
# log2.o.kim <- log2(o.kim[-1] + 1)
# log2.o.kim$gene_name <- gene.names
# log2.o.kim <- log2.o.kim %>% select(gene_name, everything())
# # check the log2-transformed data
# log2.o.kim %>% str()
# # # Save file to database
# # dbWriteTable(my.db, "ophio_kim_log2fpkm", log2.o.kim)
# 
# # zscore-log2-expression --------------------------------------------------
# # z-score the data
# gene.names <- log2.o.kim %>% pull(gene_name)
# sample.names <- names(log2.o.kim[-1])
# zscores.o.kim <- log2.o.kim %>%
#   # create a gene x exp matrix
#   select(-1) %>%
#   as.matrix() %>%
#   # use the scale function for each row to calculate z-scores
#     # scale calculates (x-mean(X))/sd(X)
#     # 1 indicates row-wise
#   apply(., 1, scale) %>%
#   # the output needs to be transposed
#   t() %>%
#   # make it a dataframe
#   as.data.frame()
# # add the column names
# names(zscores.o.kim) <- sample.names
# zscores.o.kim$gene_name <- gene.names
# zscores.o.kim <- zscores.o.kim %>% select(gene_name, everything())
# # check the z-score transformed dataset
# zscores.o.kim %>% str()
# # # Save file to database
# # dbWriteTable(my.db, "ophio_kim_zscores", zscores.o.kim)
# 
# # Save a csv with the zscores
# tbl(my.db, "ophio_kim_zscores") %>% 
#   collect() %>% 
#   na.omit() %>%
#   write.csv(.,
#             file = "./results/normalized_gene_exp/zscore/ophio_kim/ophio_kim_LD_zscores_noNAs.csv")


# 03. Ophio_kim DD --------------------------------------------------------


# o.kim <- read.csv("./results/normalized_gene_exp/raw_fpkm/ophio_kim/dd/normalized_gene_exp_ophio_kim_dd_samples.csv",
#                   header = T, stringsAsFactors = F, na.strings = c(NA, "", " "))
# # Save file to database
# dbWriteTable(my.db, "ophio_kim_DD_fpkm", o.kim)
# 
# 
# # Expressed genes ---------------------------------------------------------
# # list of all ophio_cflo "Expressed" genes (>1 FPKM during the 24h period)
# expressed <-
#   o.kim %>%
#   na.omit() %>%
#   filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>%
#   pull(gene_name) %>%
#   unique()
# # append this information to the o.kim data
# expressed.o.kim <-
#   o.kim %>%
#   select(gene_name) %>%
#   mutate(expressed = ifelse(gene_name %in% expressed,"yes","no"))
# # Save file to databse
# dbWriteTable(my.db, "ophio_kim_DD_expressed_genes", expressed.o.kim)
# 
# # log2-expression --------------------------------------------------------
# # log2-transform the data
# gene.names <- o.kim %>% pull(gene_name)
# log2.o.kim <- log2(o.kim[-1] + 1)
# log2.o.kim$gene_name <- gene.names
# log2.o.kim <- log2.o.kim %>% select(gene_name, everything())
# # check the log2-transformed data
# log2.o.kim %>% str()
# # Save file to database
# dbWriteTable(my.db, "ophio_kim_DD_log2fpkm", log2.o.kim)
# 
# # zscore-log2-expression --------------------------------------------------
# # z-score the data
# gene.names <- log2.o.kim %>% pull(gene_name)
# sample.names <- names(log2.o.kim[-1])
# zscores.o.kim <- log2.o.kim %>%
#   # create a gene x exp matrix
#   select(-1) %>%
#   as.matrix() %>%
#   # use the scale function for each row to calculate z-scores
#   # scale calculates (x-mean(X))/sd(X)
#   # 1 indicates row-wise
#   apply(., 1, scale) %>%
#   # the output needs to be transposed
#   t() %>%
#   # make it a dataframe
#   as.data.frame()
# # add the column names
# names(zscores.o.kim) <- sample.names
# zscores.o.kim$gene_name <- gene.names
# zscores.o.kim <- zscores.o.kim %>% select(gene_name, everything())
# # check the z-score transformed dataset
# zscores.o.kim %>% str()
# # Save file to database
# dbWriteTable(my.db, "ophio_kim_DD_zscores", zscores.o.kim)

# # Save a csv with the zscores
# tbl(my.db, "ophio_kim_DD_zscores") %>% 
#   collect() %>% 
#   na.omit() %>% 
#   write.csv(.,
#             file = "./results/normalized_gene_exp/zscore/ophio_kim/ophio_kim_DD_zscores_noNAs.csv",
#             row.names = F)
# 
