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


# create supplementary files ----------------------------------------------

# A1. number of 24h-rhythmic genes at GammaP < 0.05
ejtk.24 <- tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% filter(GammaP < gamma.pval) %>% collect()

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


# hierarchical clustering and heatmaps ------------------------------------

# perform hierarchical clustering of 24h-rhythmic genes into four clusters;
# plot time-course heatmaps for the clustered 24h-rhythmic geneset
# Identify the day-peaking and night-peaking clusters visually.

## Load all the rhythmic genesets 
## Note, ordered according to their p-value; highly rhythmic at the top.
# Circadian genes (period = 24h)
tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% head()
# get the gene-names for sig. 24h-rhythmic genes
rhy.24 <-
  tbl(ejtk.db, "ophio_cflo_zscores_24h") %>% 
  filter(GammaP < gamma.pval) %>% 
  select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
  select(ID) %>% pull()

# # Ultradian genes (period = 8h)
# tbl(ejtk.db, "ophio_cflo_zscores_08h") %>% head()
# ## Foragers
# rhy.8 <- 
#   tbl(ejtk.db, "ophio_cflo_zscores_08h") %>%  
#   filter(GammaP < gamma.pval) %>% 
#   select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
#   select(ID) %>% pull()
# 
# # Ultradian genes (period = 12h)
# tbl(ejtk.db, "ophio_cflo_zscores_12h") %>% head()
# ## Foragers
# rhy.12 <- 
#   tbl(ejtk.db, "ophio_cflo_zscores_12h") %>% 
#   filter(GammaP < gamma.pval) %>% 
#   select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
#   select(ID) %>% pull()

## load zscore dataset
zscore.24h <- data.db %>% tbl(., "ophio_cflo_zscores") %>% collect()

# Filter the zscores to keep only circadian genes
zscore.24h <-
  zscore.24h %>% 
  filter(gene_name %in% rhy.24) %>% 
  as.data.frame()

# Set genes as rownames and convert it into a matrix
rownames(zscore.24h) = zscore.24h$gene_name
zscore.24h <- as.matrix(zscore.24h[-1])


# Hierarchical clustering of the for24 and nur24 genesets
my_hclust_gene <- hclust(dist(zscore.24h), method = "complete")


# Make annotations for the heatmaps
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4) # four clusters
my_gene_col <- data.frame(cluster = my_gene_col)


# I’ll add some column annotations and create the heatmap.
# Annotations for:
# 1. Is the sample collected during the light or dark phase? 
my_sample_col <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col) <- colnames(zscore.24h)

# Manual color palette
my_colour = list(
  phase = c(light = "#F2E205", dark = "#010440"),
  cluster = viridis::cividis(100)[c(10,90,60,30)])

# Color scale
my.breaks = seq(-3, max(zscore.24h), by=0.06)

# Let's plot!
# rhy.heat <- 
  pheatmap(zscore.24h, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col, 
                         annotation_col = my_sample_col,
                         cutree_rows = 4,
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



# run enrichment for day/night-peaking clusters ---------------------------

# Run enrichment for diurnal (24h-rhythmic) genes
# that peak during the day (day-peaking clusters) and night (night-peaking clusters);
  # save the results as csvs & 
  # plot the enrichment results.

# from eye-balling the heatmap from above
  # day-peaking cluster: cluster-3
  # night-peaking cluster: cluster-1
  
## day-peaking | cluster 3 ##
rhy.24.daypeaking.cluster3 <-
  my_gene_col %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster == 3) %>%
  pull(gene) %>%
  # run enrichment analysis
  go_enrichment(., 
                org = "ophio_cflo", 
                bg = "expressed") # enrichment against all expressed ophio_cflo genes
# view the results
rhy.24.daypeaking.cluster3 %>% view()

## night-peaking | cluster 1 ##
rhy.24.nightpeaking.cluster1 <- 
  my_gene_col %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster == 1) %>%
  pull(gene) %>%
  go_enrichment(., 
                org = "ophio_cflo", 
                bg = "expressed")
# view the results
rhy.24.nightpeaking.cluster1 %>% view()
