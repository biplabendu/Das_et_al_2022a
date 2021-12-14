set.seed(420)
rm(list = ls())
#
# set working dirictory
# setwd("~/Dropbox/Ant-fungus/02_git/Git_Das_folder2/Das_et_al_2022a")
#
## Load packages ----------
pacman::p_load(pheatmap, dendextend, tidyverse, viridis, ggthemes)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
pacman::p_load(glue)
#
# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
#
#' set path to your working directory
path_to_repo = "/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a"

# script name
script.name = "03_compare_circadianGCN_ophio-cflo_v_ophio-kim"

species <- 'ophio_cflo'

#' load functions
# customized theme for publication quality figures
source(paste0(path_to_repo,"/functions/theme_publication.R"))
# function to perform enrichment analysis
source(paste0(path_to_repo,"/functions/enrichment_analysis.R"))
# function to plot z-scores (Cflo genes only)
source(paste0(path_to_repo,"/functions/plot_zscores.R"))

path_to_repo = "/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a"

# Connect to db's
rhy.db <- dbConnect(RSQLite::SQLite(),
                    "./data/databases/TC6_fungal_ejtk.db")


## Overview/Goals
 # - check homologs of night peaking genes in ophio_cflo and ophio_kim

## Step 1: Identify the homologous genes

### 1.1 Load, clean and merge data

# Dataset: Ophio cflo data (as csv file), containing information about rhythmic genes and ophio kim homologs. Ophio kim data (stored in DB), containing information about rhythmic genes in Ophio_kim


# load cflo data with the with the following info: Robin gene ID's of cflo, homologs of kim in cflo, 24h rhythmic genes in cflo and the gammaP value of the 24h rhythmic genes 
cflo.data <- read.csv(glue('{path_to_repo}/data/ophio_cflo_TC6_data.csv')) %>% 
  select(ophio_cflo_ID = gene_ID_robin, ophio_kim_ID = ophio_kim_homolog, rhythmic_24h_ophio_cflo = rhythmic_24h)

# load the ophio kim data with the following info: Robin's gene ID's, 24h rhythmic genes in ophio_kim and their GammaP value
kim.data <- rhy.db %>% 
  tbl('ophio_kim_LD_rhythmic_genes_24h') %>% 
  select(rhythmic_24h_ophio_kim = rhythmic, ophio_kim_ID = gene_ID_robin) %>% 
  collect()


#### Other try

# unique one-to-one homologous okim genes
orthos.kim <-
  cflo.data %>% 
  select(-rhythmic_24h_ophio_cflo) %>% 
  distinct() %>% 
  na.omit() %>% 
  group_by(ophio_kim_ID) %>% 
  summarise(num = n()) %>% 
  filter(num == 1) %>%
  pull(ophio_kim_ID) %>% 
  as.character()


# Merge the data on othologs of ophio kim with ophio cflo
# data <- left_join(kim.data, orthos, by='ophio_kim_ID')
data <-
  cflo.data %>% 
  filter(ophio_kim_ID %in% orthos.kim) %>%
  left_join(kim.data, by="ophio_kim_ID")


# write.csv(data, 'unique_orthos_cflo_kim.csv')

############
# Night-peaking ophio.kim (134 genes)
bar <- 
  read_csv(glue('{path_to_repo}/data/ophio_kim_2017.csv')) %>% 
  filter(Clustered_Activity == 'night-active')  %>% 
  unique() %>% 
  left_join(data, by='ophio_kim_ID') 

write_csv(bar, 'foo.csv')

%>% 
  pull(ophio_cflo_ID)

length(bar)

# Night-peaking ophio.cflo (833 genes) 
foo <- 
  read_csv(glue('{path_to_repo}/data/rdata/ophio_cflo_hierachical_clusters.csv')) %>% 
  filter(cluster == '1') %>% 
  pull(gene_ID_robin)

length(foo)

# rhytmic genes in ophio cflo and kim as ophio cflo ID's (342)
both_rhy4 <- data %>% 
  filter( rhythmic_24h_ophio_kim == 'yes' & rhythmic_24h_ophio_cflo == 'yes') %>% 
  unique() %>% 
  pull(ophio_cflo_ID)

# non rhythmic genes in ophio kim and cflo, as ophio cflo ID's (3506)
not_rhy4 <- data %>% 
  filter( rhythmic_24h_ophio_kim == 'no' & rhythmic_24h_ophio_cflo == 'no')%>% 
  unique() %>% 
  pull(ophio_cflo_ID)


#################
# overlap between rhy_cflo and rhy_kim (342)
x1 <- intersect(foo, bar) %>% length()
# genes in cflo_rhy but not in kim_rhy (1762)
x2 <- setdiff(foo, bar) %>% length()
# genes in kim_rhy and not in cflo_rhy (683)
x3 <- setdiff(bar, foo) %>% length()
# genes not in a and not in b, or just total amount of genes? Just picked 7455
x4 <- 7455

# make the table
con.table <- data.frame(c1 = c(x1,x2),
                        c2 = c(x3,x4))

fisher.test(con.table)


#################
# Overlap DEG signalP and rhytmic genes ophio
# overlap between rhy_cflo and rhy_kim (342)
x1 <- 124
# genes in cflo_rhy but not in kim_rhy (1762)
x2 <- 320
# genes in kim_rhy and not in cflo_rhy (683)
x3 <- 2285
# genes not in a and not in b, or just total amount of genes? Just picked 7455
x4 <- 2481

# make the table
con.table <- data.frame(c1 = c(x1,x2),
                        c2 = c(x3,x4))

fisher.test(con.table)
