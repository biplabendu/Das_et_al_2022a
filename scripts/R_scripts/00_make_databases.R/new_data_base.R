########## TC6 Fungal Database ############
#
# This script is written to create a databse containing all the data of a 24h time-course experiment
#
# Housekeeping ---------------------------------------------------------------
#
set.seed(420)
rm(list = ls())
#
## Load packages --
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
pacman::p_load(glue)
#
# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
#
## set parameters and thresholds
#
# Specify your path
path <- "/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a/"
#
# Specify species
species <- 'ophio_cflo'

# 0. Create the database / connect to the database ----------------------------------
my.db <- dbConnect(RSQLite::SQLite(),
                   glue("{path}/data/databases/new_TC6_fungal_data.db"))


# 1. FPKM table  ----------------------------------------------------------
#
# Load in the expression data
exp.val <- read.csv(glue("{path}/results/normalized_gene_exp/raw_fpkm/{species}/normalized_gene_exp_{species}_all_samples.csv"),
                 header = T, stringsAsFactors = F, na.strings = c(NA, "", " ")) %>%
  as.data.frame()

# also make a colum for the gene ID's from Robin
#
# Load in the gene IDs
gene.IDs <- read.csv(glue("{path}/data/input/{species}/{species}_gene_names_robin_ncbi.csv"),
                     header = T, stringsAsFactors = F, na.strings = c(NA, "", " ")) %>%
  as.data.frame() 

# Join the data
data <- left_join(exp.val, gene.IDs, by = 'gene_ID_ncbi')
# rearange cols
data <- data[,c(16,1:15)]

# Save file to database
dbWriteTable(my.db, glue("{species}_fpkm"), data)

# check if table is added in the database (by checking all tables)
src_dbi(my.db)

# A. Expressed genes ---------------------------------------------------------

# list of all beau "Expressed" genes (>1 FPKM during the 24h period)
expressed <-
  my.db %>%
  tbl(glue("{species}_fpkm")) %>% 
  select(-end,
         -start)%>% 
  collect() %>%
  na.omit() %>% # Remove NA's
  filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>% # expression > 1 FPKM 
  pull(gene_ID_ncbi) %>% # get the gene names of the expressed genes
  unique() # check wheter there are duplicates

# append this information to the beau data
expressed.sp <-
  data %>%
  select(gene_ID_ncbi, gene_ID_robin) %>%
  mutate(expressed = ifelse(gene_ID_ncbi %in% expressed,"yes","no"))

# Save file to databse
dbWriteTable(my.db, glue("{species}_expressed_genes"), expressed.sp)

# show tables in DB
src_dbi(my.db)

# # B. log2-expression --------------------------------------------------------
# # log2-transform the data
# gene.names <- beau %>% pull(gene_name)
# log2.beau <- log2(beau[-1] + 1)
# log2.beau$gene_name <- gene.names
# log2.beau <- log2.beau %>% select(gene_name, everything())
# # check the log2-transformed data
# log2.beau %>% str()
# # # Save file to database
# dbWriteTable(my.db, "beau_log2fpkm", log2.beau)
# 
# # C. score-log2-expression --------------------------------------------------
# # z-score the data
# gene.names <- log2.beau %>% pull(gene_name)
# sample.names <- names(log2.beau[-1])
# zscores.beau <- log2.beau %>%
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
# names(zscores.beau) <- sample.names
# zscores.beau$gene_name <- gene.names
# zscores.beau <- zscores.beau %>% select(gene_name, everything())
# # check the z-score transformed dataset
# zscores.beau %>% str()
# # Save file to database
# dbWriteTable(my.db, "beau_zscores", zscores.beau)
# 
# # Save a csv with the zscores
# beau.zscore <- tbl(my.db, "beau_zscores") %>% collect()
# 
# beau.zscore.noNAs <-
#   beau.zscore %>%
#   na.omit()
# 
# write.csv(beau.zscore.noNAs,
#           file = "./results/normalized_gene_exp/zscore/beau/beau_zscores_noNas.csv",
#           row.names = F)
