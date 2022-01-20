########## TC6 Fungal Database ############
#
# This script is written to create a databse containing all the data of a obtained during Time Course experiment 6
# This database can constantly be updated to add more columns to the excisting tables
# The end goal is to have one big table per species, which contains all the data, that can be easily written to a cvs and used for sharing
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
src_dbi(my.db)


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

# B. 
