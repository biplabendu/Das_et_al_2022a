### Housekeeping
set.seed(420)
rm(list = ls())

#' Load the libraries
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted, WGCNA, igraph)
pacman::p_load(patchwork)
pacman::p_load(glue)

#' set conflict preference
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
conflict_prefer("hclust", "flashClust")
conflict_prefer("simplify", "igraph")

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


# Connect to db's
rhy.db <- dbConnect(RSQLite::SQLite(),
                    "./data/databases/TC6_fungal_ejtk.db")

### DATA load

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

# ophio cflo rhythmic genes as ophio cflo gene ID's (2104 genes)
cflo_rhy24h <- data %>% 
  filter(rhythmic_24h_ophio_cflo == 'yes') %>% 
  unique() %>% 
  pull(ophio_cflo_ID)

# ophio cflo rhythmic genes as ophio gene cflo ID's (1025)
kim_rhy24h <- data %>% 
  filter(rhythmic_24h_ophio_kim == 'yes') %>% 
  unique() %>% 
  pull(ophio_cflo_ID)

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



### Constituently table


# overlap between rhy_cflo and rhy_kim (342)
x1 <- intersect(cflo_rhy24h,kim_rhy24h) %>% length()
# genes in cflo_rhy but not in kim_rhy (1762)
x2 <- setdiff(cflo_rhy24h, kim_rhy24h) %>% length()
# genes in kim_rhy and not in cflo_rhy (683)
x3 <- setdiff(kim_rhy24h, cflo_rhy24h) %>% length()
# genes not in a and not in b, or just total amount of genes? Just picked 7455
x4 <- 7455

# make the table
con.table <- data.frame(c1 = c(x1,x2),
                        c2 = c(x3,x4))