# Housekeeping ---------------------------------------------------------------
#
set.seed(420)
rm(list = ls())
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
source("./functions/enrichment_analysis.R")
#
## set parameters and thresholds --------
#
# gamma-pvalue threshold for inferring rhythmicity
gamma.pval = 0.05
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


# Specify the focal sample for analysis -----------------------------------

sample.name = "beau"

# 01. General patterns of gene expression ---------------------------------

# A1: genes that have NO expression (FPKM == 0 at all time points)
not.expressed <-
  tbl(data.db, paste0(sample.name ,"_fpkm")) %>% 
  collect() %>% 
  filter_at(vars(starts_with("Z")), all_vars(. == 0)) %>%
  pull(gene_name)

# Write all the non-expressed genes to a file

write(not.expressed, file = paste0('./results/',{sample.name},'_not_expressed_list.txt'),
     sep = " ")

# A2: run enrichment (make plot of enrichment found of non-expressed genes)
not.expressed %>% 
  go_enrichment(., 
                org = "beau", 
                bg = 'all') %>%  # enrichment against all ophio_cflo genes in the genome
  go_enrichment_plot(clean = "no")
  
# B: genes that are expressed (FPKM > 1 for at least one time point)
expressed <- 
  tbl(data.db, paste0(sample.name,"_expressed_genes")) %>% 
  filter(expressed=="yes") %>% 
  collect() %>% 
  pull(gene_name)


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
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2) # k=  clusters
my_gene_col <- data.frame(cluster = my_gene_col)


# I’ll add some column annotations and create the heatmap.
# Annotations for:
# 1. Is the sample collected during the light or dark phase? 
my_sample_col <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col) <- colnames(zscore.rhy)

# Manual color palette
my_colour = list(
  phase = c(light = "#F2E205", dark = "#010440"),
  cluster = viridis::cividis(100)[c(40,60)])

# Color scale
my.breaks = seq(min(zscore.rhy), max(zscore.rhy), by=0.1)
# my.breaks = seq(min(zscore.rhy), max(zscore.rhy), by=0.06)

# Let's plot!
rhy.heat <-
  pheatmap(zscore.rhy, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col, 
                         annotation_col = my_sample_col,
                         cutree_rows = 2, # OG was 4
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
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
name.path.file <- paste0('./results/heatmap/rhy_heatmap_',sample.name,'_',period,'h_2clusters.pdf')
save_pheatmap_pdf(rhy.heat, name.path.file)
  
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
                org = "beau",
                bg = expressed) # enrichment against all expressed ophio_cflo genes
# view the results
rhy.24.daypeaking.cluster3 %>% view()

## night-peaking | cluster 1 ##
rhy.24.nightpeaking.cluster1 <-
  my_gene_col %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster == 1) %>%
  pull(gene) %>%
  go_enrichment(.,
                org = "beau",
                bg = "expressed")
# view the results
rhy.24.nightpeaking.cluster1 %>% view()

