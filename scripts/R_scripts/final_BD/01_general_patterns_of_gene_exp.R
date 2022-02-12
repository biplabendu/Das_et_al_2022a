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
pacman::p_load(glue)
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
# set path to your working directory
path_to_repo = "/Users/biplabendudas/Documents/GitHub/Das_et_al_2022a"
# gamma-pvalue threshold for inferring rhythmicity
gamma.pval = 0.05
#
#
## Path to save files
# # for supplementary files
# supp.path="~/University\ of\ Central\ Florida/Charissa\ De\ Bekker\ -\ Ant-Fungus-Clock-Interactions/04_manuscript/03_supplementary_files/"



# 00. Load Databases -----------------------------------------------------------
#
# 1. TC6_ejtk.db
# Desc: This database contains all ejtk-output for TC6
ejtk.db <- dbConnect(RSQLite::SQLite(),
                     "./data/databases/TC6_fungal_ejtk.db")
# which tables are in the database
src_dbi(ejtk.db)
#
# 2. TC6_data.db
data.db <- dbConnect(RSQLite::SQLite(),
                     "./data/databases/TC6_fungal_data.db")
src_dbi(data.db)
#
##

# 01. General patterns of gene expression ---------------------------------
#
# specify sample for analysis
sample.name <- c("beau", "ophio_cflo", "ophio_kim")
#
# number of all genes
all.genes <- list()
for (i in 1:length(sample.name)) {
  all.genes[[i]] <- tbl(data.db, paste0(sample.name[[i]] ,"_fpkm")) %>%  
    collect()
  
  writeLines(paste("Number of genes in", sample.name[[i]], ":", nrow(all.genes[[i]])))
}

# A1: genes that have NO expression (FPKM == 0 at all time points)
not.expressed <- list()
for (i in 1:length(sample.name)) {
  not.expressed[[i]] <-
    tbl(data.db, paste0(sample.name[[i]] ,"_fpkm")) %>% 
    collect() %>% 
    filter_at(vars(starts_with("Z")), all_vars(. == 0)) %>%
    pull(gene_name)
  
  # How many genes are not expressed?
  writeLines(paste("n(genes-NOT-EXPRESSED) in", sample.name[[i]], ":", length(not.expressed[[i]])))
  
}

# A2: run enrichment (make plot of enrichment found of non-expressed genes)
for (i in 1:length(sample.name)) {
  writeLines(paste("running GO enrichment for NOT-EXPRESSED genes in", sample.name[[i]]))
  # run enrichment
  not.expressed[[i]] %>% 
    go_enrichment(., 
                  org = sample.name[[i]], 
                  bg = 'all') %>% # enrichment against all ophio_cflo genes in the genome
    
    # # pull gene names for a given GO term
    # separate_rows(., gene_name, sep = ", ") %>%
    # filter(GO == "GO:0009405") %>% # pathogenesis
    # # filter(GO == "GO:0090729") %>% # toxin activity
    # # filter(GO == "GO:0044419") %>% # interspecies interaction between organisms
    # # filter(GO == "GO:0020037") %>% # heme binding
    # pull()
    
    go_enrichment_plot(clean = "no") %>% 
    print()
  
}

# B: genes that are expressed (FPKM > 1 for at least one time point)
expressed <- list()
for (i in 1:length(sample.name)) {
  expressed[[i]] <- 
    tbl(data.db, paste0(sample.name[[i]],"_expressed_genes")) %>% 
    filter(expressed=="yes") %>% 
    collect() %>% 
    pull(gene_name) 
  
  # How many genes are expressed?
  writeLines(paste("n(EXPRESSED) in", sample.name[[i]], ":", length(expressed[[i]])))
}



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

## 
rhy <- list()
for (i in 1:2) {
  rhy[[i]] <-
    tbl(ejtk.db, paste0(sample.name[[i]],"_zscores_",period,'h')) %>%
    filter(GammaP < gamma.pval) %>%
    select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
    select(ID) %>% pull()
  
  # How many genes are rythmic?
  writeLines(paste0("n(rhythmic-",period, "h) in ", sample.name[[i]], " : ", length(rhy[[i]])))
}

## initialise lists to hold input and output of the hierarchical clustering
zscore.dat <- list() # zscore data (input)
my_gene_col <- list() # cluster identity for each rhythmic gene (output)
rhy.heat <- list() # pheatmap that can be saved/plotted (output)
## run clustering and plot
for (i in 1:2) {
  ## load zscore dataset
  zscore.dat[[i]] <- data.db %>% tbl(., paste0(sample.name[[i]],"_zscores")) %>% collect()
  
  # Filter the zscores to keep only rhythmic genes
  zscore.rhy <-
    zscore.dat[[i]] %>% 
    filter(gene_name %in% rhy[[i]]) %>% 
    as.data.frame()
  
  # Set genes as rownames and convert it into a matrix
  rownames(zscore.rhy) = zscore.rhy$gene_name
  zscore.rhy <- as.matrix(zscore.rhy[-1])
  
  
  # Hierarchical clustering of the genesets
  my_hclust_gene <- hclust(dist(zscore.rhy), method = "complete")
  
  
  # Make annotations for the heatmaps
  my_clusters <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4) # k=  clusters
  my_gene_col[[i]] <- data.frame(cluster = my_clusters)
  
  
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
  rhy.heat[[i]] <-
    pheatmap(zscore.rhy, show_rownames = F, show_colnames = F,
             annotation_row = my_gene_col[[i]], 
             annotation_col = my_sample_col,
             cutree_rows = 4, # OG was 4
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
  
}


# Beau-24h-clusters -------------------------------------------------------

# Cluster 1
writeLines("Cluster 1; n = 767 genes")
writeLines("a ton of enriched GO terms, showing terms at FDR = 0.1%")
my_gene_col[[1]] %>% 
  rownames_to_column("gene") %>%
  filter(cluster %in% c("1")) %>%
  pull(gene) %>%
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[1]],
                bg = expressed[[1]]) %>%
  
  # separate_rows(., gene_name) %>%
  
  go_enrichment_plot(clean = "no", fdr = 0.1)

# Cluster 2
writeLines("Cluster 2; n = 550 genes")
my_gene_col[[1]] %>% 
  rownames_to_column("gene") %>%
  filter(cluster %in% c("2")) %>% 
  pull(gene) %>%
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[1]],
                bg = expressed[[1]]) %>%
  
  # separate_rows(., gene_name) %>% view()
  
  go_enrichment_plot(clean = "no")

# # Cluster 3
writeLines("Cluster 3; n = 337 genes; no enriched GO terms")
# my_gene_col[[1]] %>% 
#   rownames_to_column("gene") %>%
#   filter(cluster %in% c("3")) %>%
#   pull(gene) %>% 
#   
#   # stacked.zplot_tc7() %>%
#   # keep(names(.) %in% c("control","ophio")) %>%
#   # multi.plot(rows = 2,cols = 1)
#   
#   go_enrichment(.,
#                 function.dir = path_to_repo,
#                 org = sample.name[[1]],
#                 bg = (expressed[[1]])) %>%
#   
#   # separate_rows(., gene_name) %>% view()
#   
#   go_enrichment_plot(clean = "no")

# Cluster 4
writeLines("Cluster 4; n = 218 genes")
my_gene_col[[1]] %>% 
  rownames_to_column("gene") %>%
  filter(cluster %in% c("4")) %>%
  pull(gene) %>% 
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[1]],
                bg = (expressed[[1]])) %>%
  
  # separate_rows(., gene_name) %>% view()
  
  go_enrichment_plot(clean = "no")


# Ocflo-24h-clusters ------------------------------------------------------

# Cluster 1
writeLines("Cluster 1; n = 833 genes")
writeLines("a ton of enriched GO terms, showing terms at FDR = 0.1%")
png(paste0(path_to_repo, "/results/figures/BD/", sample.name[2], "/", sample.name[2],"_rhy24_cluster1.png"),
    width = 20, height = 35, units = "cm", res = 300)
my_gene_col[[2]] %>% 
  rownames_to_column("gene") %>%
  filter(cluster %in% c("1")) %>%
  pull(gene) %>%
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[2]],
                bg = expressed[[2]]) %>% 
  
  # separate_rows(., gene_name) %>% 
  
  go_enrichment_plot(clean = "no", fdr = 0.1)
dev.off()

# Cluster 2
writeLines("Cluster 2; n =  genes")
my_gene_col[[2]] %>% 
  rownames_to_column("gene") %>%
  filter(cluster %in% c("2")) %>% 
  pull(gene) %>%
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[2]],
                bg = expressed[[2]]) %>%
  
  # separate_rows(., gene_name) %>% view()
  
  go_enrichment_plot(clean = "no")

# # Cluster 3
writeLines("Cluster 3; n = 354 genes")
my_gene_col[[2]] %>%
  rownames_to_column("gene") %>%
  filter(cluster %in% c("3")) %>%
  pull(gene) %>%
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[2]],
                bg = (expressed[[2]])) %>%
  
  # separate_rows(., gene_name) %>% view()
  
  go_enrichment_plot(clean = "no")

# Cluster 4
writeLines("Cluster 4; n = 633 genes")
my_gene_col[[2]] %>% 
  rownames_to_column("gene") %>%
  filter(cluster %in% c("4")) %>%
  pull(gene) %>% 
  
  # stacked.zplot_tc7() %>%
  # keep(names(.) %in% c("control","ophio")) %>%
  # multi.plot(rows = 2,cols = 1)
  
  go_enrichment(.,
                function.dir = path_to_repo,
                org = sample.name[[2]],
                bg = (expressed[[2]])) %>%
  
  # separate_rows(., gene_name) %>% view()
  
  go_enrichment_plot(clean = "no")


