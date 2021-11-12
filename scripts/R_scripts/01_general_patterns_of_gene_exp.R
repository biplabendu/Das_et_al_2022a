# Housekeeping ---------------------------------------------------------------
#
set.seed(420)
rm(list = ls())
#
# set working dirictory
# setwd("~/Dropbox/Ant-fungus/02_git/Git_Das_folder2/Das_et_al_2022a")
#
## Load packages ----------
pacman::p_load(pheatmap, dendextend, tidyverse, viridis, ggthemes)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
#
# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
#
## load functions ---------
# functions for enrichment
source("./functions/enrichment_analysis.R")
# customized theme for publication quality figures
source("../Das_et_al_2022b/functions/theme_publication.R")
#
## set parameters and thresholds --------
#
# gamma-pvalue threshold for inferring rhythmicity
gamma.pval = 0.05
#
#
## Path to save files
# # for supplementary files
# supp.path="~/University\ of\ Central\ Florida/Charissa\ De\ Bekker\ -\ Ant-Fungus-Clock-Interactions/04_manuscript/03_supplementary_files/"



# 00. Load Databases -----------------------------------------------------------
#
# 1. TC7_ejtk.db
# Desc: This database contains all ejtk-output for TC7
ejtk.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC6_fungal_ejtk.db")
# which tables are in the database
src_dbi(ejtk.db)
#
# 2. TC7_data.db
data.db <- dbConnect(RSQLite::SQLite(),
                     "./data/databases/TC6_fungal_data.db")
src_dbi(data.db)
#
##

# 01. General patterns of gene expression ---------------------------------
#
# specify sample for analysis
sample.name <- 'ophio_cflo'
#
# number of all genes
all.genes <- tbl(data.db, paste0(sample.name ,"_fpkm")) %>%  
  collect() %>% 
  pull(gene_name)

length(all.genes)

# A1: genes that have NO expression (FPKM == 0 at all time points)
not.expressed <-
  tbl(data.db, paste0(sample.name ,"_fpkm")) %>% 
  collect() %>% 
  filter_at(vars(starts_with("Z")), all_vars(. == 0)) %>%
  pull(gene_name)

# How many genes are not expressed?
length(not.expressed)

#
# Write all the non-expressed genes to a file
# write(not.expressed, file = paste0('./results/',{sample.name},'_not_expressed_list.txt'), sep = " ")

# A: run enrichment (make plot of enrichment found of non-expressed genes)
not.expressed %>% 
  go_enrichment(., 
                org = sample.name, 
                bg = 'all') %>%  # enrichment against all ophio_cflo genes in the genome
  go_enrichment_plot(clean = "no")
  
# B: genes that are expressed (FPKM > 1 for at least one time point)
expressed <- 
  tbl(data.db, paste0(sample.name,"_expressed_genes")) %>% 
  filter(expressed=="yes") %>% 
  collect() %>% 
  pull(gene_name)
#
# How many genes are expressed?
length(expressed)



# 02. Diurnal rhythms in gene expression ----------------------------------


# create supplementary files ----------------------------------------------

# A1. number of 24h-rhythmic genes at GammaP < 0.05
# ejtk.24 <- tbl(ejtk.db, paste0(sample.name,"_zscores_24h")) %>% filter(GammaP < gamma.pval) %>% collect()

# A2. save the ejtk-results (ALL GENES) as a supplementary file
# tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% 
#   collect() %>% 
#   write.csv(.,
#             file = paste0(supp.path,"S1_ejtk_24h_ophio_cflo_all_genes.csv"),
#             row.names = F)

# A3. save the ejtk-results (rhy24 genes only) as a supplementary file
# tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% 
#   filter(GammaP < gamma.pval) %>% 
#   collect() %>% 
#   write.csv(.,
#             file = paste0(supp.path,"S2_ejtk_24h_ophio_cflo_rhy_genes.csv"),
#             row.names = F)


# 03. Hierarchical clustering and heatmaps ------------------------------------

# perform hierarchical clustering of 24h-rhythmic genes into four clusters;
# plot time-course heatmaps for the clustered 24h-rhythmic geneset
# Identify the day-peaking and night-peaking clusters visually.

## Load all the rhythmic genesets 
## Note, ordered according to their p-value; highly rhythmic at the top.
#
# Choose period
period = '24'

# Ultradian genes (period = 8h)
## 
rhy <-
  tbl(ejtk.db, paste0(sample.name,"_zscores_",period,'h')) %>%
  filter(GammaP < gamma.pval) %>%
  select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
  select(ID) %>% pull()

# How many genes are rythmic?
length(rhy)
print(paste0("genes are rythmic expression during", period,"hours"))

## load zscore dataset
zscore.dat <- data.db %>% tbl(., paste0(sample.name,"_zscores")) %>% collect()

# Filter the zscores to keep only rhythmic genes
zscore.rhy <-
  zscore.dat %>% 
  filter(gene_name %in% rhy) %>% 
  as.data.frame()


# Set genes as rownames and convert it into a matrix
rownames(zscore.rhy) = zscore.rhy$gene_name
zscore.rhy <- as.matrix(zscore.rhy[-1])


# Hierarchical clustering of the genesets
my_hclust_gene <- hclust(dist(zscore.rhy), method = "complete")


# Make annotations for the heatmaps
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4) # k=  clusters
my_gene_col <- data.frame(cluster = my_gene_col)


# I’ll add some column annotations and create the heatmap.
# Annotations for:
# 1. Is the sample collected during the light or dark phase? 
my_sample_col <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col) <- colnames(zscore.rhy)

# Manual color palette
my_colour = list(
  phase = c(light = "#F2E205", dark = "#010440"),
  cluster = viridis::cividis(100)[c(10,90,60,30)])

# Color scale
my.breaks = seq(min(zscore.rhy), max(zscore.rhy), by=0.1)
# my.breaks = seq(min(zscore.rhy), max(zscore.rhy), by=0.06)

# Let's plot!
rhy.heat <-
  pheatmap(zscore.rhy, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col, 
                         annotation_col = my_sample_col,
                         cutree_rows = 1, # OG was 4
                         cutree_cols = 2,
                         annotation_colors = my_colour,
                         border_color=FALSE,
                         cluster_cols = F,
                         breaks = my.breaks,
                         ## color scheme borrowed from: 
                         color = inferno(length(my.breaks) - 1),
                         # treeheight_row = 0, 
                         # treeheight_col = 0,
                         # remove the color scale or not
                         # main = paste0("Foragers - circadian genes \n (n=", nrow(cflo.rhy.exp.for), " genes)"),
                         ## annotation legend
                         annotation_legend = T,
                         ## Color scale
                         legend = T)

# To save the heatmap to a pdf, run this code. For this to work make sure the heatmap is stored in the variable rhy.heat
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# name.path.file <- paste0('./results/heatmap/rhy_heatmap_',sample.name,'_',period,'h_2clusters.pdf')
# save_pheatmap_pdf(rhy.heat, name.path.file)
  
# run enrichment for day/night-peaking clusters ---------------------------

# Run enrichment for diurnal (24h-rhythmic) genes
# that peak during the day (day-peaking clusters) and night (night-peaking clusters);
  # save the results as csvs &
  # plot the enrichment results.

# Which clusters are day and night peaking? 
# from eye-balling the heatmap from above
  # day-peaking cluster: cluster-3
  # night-peaking cluster: cluster-1
# 
# ## day-peaking | cluster 3 ##
# rhy.daypeaking.cluster <-
#   my_gene_col %>%
#   rownames_to_column(var = "gene") %>%
#   filter(cluster %in% c(2,4)) %>% # %in% c=(1,2) for filter in 2, == for one # NAME here the cluster
#   pull(gene) #%>%

# # save day peaking cluster to a file
# cluster.file.name.D <- paste0('./results/temp_files/',sample.name,'_',period,'h.txt')
# lapply(rhy.daypeaking.cluster, write, cluster.file.name.D, append=TRUE, ncolumns=1000)
# length(rhy.daypeaking.cluster)

# # run enrichment analysis
#   go_enrichment(rhy.daypeaking.cluster,
#                 org = sample.name,
#                 bg = expressed) # enrichment against all expressed ophio_cflo genes
# # view the results
# rhy.daypeaking.cluster %>% view()
# 
# # Plotting the enriched GOs for day-peaking clusters
# rhy.daypeaking.cluster %>%
#   go_enrichment_plot(clean = "no",
#                      # function.dir = path_to_repo,
#                      fdr = 5)
# # # view the results
# # rhy.24.daypeaking.cluster.3.4 %>% head() %>% view()
# 
# 
# ## night-peaking | cluster 1 ##
# rhy.nightpeaking.cluster <-
#   my_gene_col %>%
#   rownames_to_column(var = "gene") %>%
#   filter(cluster == 1) %>%
#   pull(gene)

# %>%
#   go_enrichment(.,
#                 org = "beau",
#                 bg = expressed)
# 
# # view the results and save GO-terms enriched with their p-value to csv
# rhy.nightpeaking.cluster %>% 
#   filter(adj_pVal < 0.05) %>% 
#   filter(over_under == "over") %>% 
#   select(GO, adj_pVal) %>% 
#   write.csv(., file = "./results/billu_beau_night_peaking.csv", row.names = F)
#   # view()
# #
# # plotting the enriched GOs for night-peaking clusters")
# # plot the enriched GOs
# rhy.nightpeaking.cluster %>%
#   go_enrichment_plot(clean = "no", 
#                      # function.dir = path_to_repo,
#                      fdr = 1)

################ Export for enrichment on robins site -------------------------
#get robin names
names_robin_ncbi <- read_csv("data/input/ophio_cflo/ophio_cflo_gene_names_robin_ncbi.csv")

merge <- 
  my_gene_col %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster == 2) # %in% c(1,3)) # %in% c=(1,2) for filter in 2, == for one  

hoi <- left_join(merge, names_robin_ncbi, by=c('gene' = 'gene_ID_ncbi'))
cluster.file.name <- paste0('./results/temp_files/',sample.name,'_',period,'_cluster_Ye.csv')
write.csv(hoi, cluster.file.name)

# cluster.file.name.N <- paste0('./results/peaking_clusters/night_peaking_Ncbi_IDs_',sample.name,'_',period,'h.txt')
# lapply(rhy.nightpeaking.cluster, write, append=TRUE, cluster.file.name.N, ncolumns=2000)
# length(rhy.nightpeaking.cluster)

####### File with all gene names and their cluster from hierachical clustering

new <- 
  my_gene_col %>%
  rownames_to_column(var = "gene") 

## add gene_ID_robin

#first get rid of extra cols
names_robin_ncbi <- read_csv("data/input/ophio_cflo/ophio_cflo_gene_names_robin_ncbi.csv")

names_robin_ncbi$start=NULL
names_robin_ncbi$end=NULL

hoi2 <- left_join(names_robin_ncbi, new, by=c('gene_ID_ncbi'= 'gene'))
file.name <- paste0('./results/temp_files/',sample.name,'_',period,'_hierachical_clusters.csv')
write.csv(hoi2, file.name)

###################################### ---------------------------
