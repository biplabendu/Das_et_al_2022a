set.seed(420)
rm(list = ls())

#' Load the libraries
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted, igraph)
pacman::p_load(patchwork, glue)

#' set conflict preference
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
conflict_prefer("hclust", "flashClust")

#' set path to your working directory
path_to_repo = "."

#' load functions
# customized theme for publication quality figures
source(paste0(path_to_repo,"/functions/theme_publication.R"))
# function to perform enrichment analysis
source(paste0(path_to_repo,"/functions/enrichment_analysis.R"))
# function to plot z-scores (Cflo genes only)
source(paste0(path_to_repo,"/functions/plot_zscores.R"))


# loading database which contains data for Das and de Bekker 2021 (bioRxiv)
db <- dbConnect(RSQLite::SQLite(), paste0(path_to_repo,"/data/databases/new_TC6_fungal_data.db"))
src_dbi(db)

# specify sample name
sample.name <- "ophio_cflo"
file <- read_csv(glue('./data/{sample.name}_TC6_data.csv'))

# extract the (gene-expr X time-point) data
dat <-
  db %>%
  tbl(., paste0(sample.name ,"_fpkm")) %>%
  select(gene_name = gene_ID_robin, everything()) %>%
  select(-c(start,end, gene_ID_ncbi)) %>% 
  collect()


# Which genes are expressed throughout the day in forager heads?
# count the number of time points that has ≥ 1 FPKM
n.expressed <- apply(dat[-1], 1, function(x) sum(x >= 1))
# subset the data and only keep the filtered genes
dat <- dat[which(n.expressed >=6),]

datExpr = as.data.frame(t(log2(dat[-c(1)]+1)))
names(datExpr) = dat$gene_name
rownames(datExpr) = names(dat)[-c(1)]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#######


pacman::p_load(GeneOverlap)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf
mergedMEs=read_csv('./results/networks/ophio_cflo_gene_IDs_and_module_identity.csv')

colors <- mergedMEs$module_identity %>% unique()
module_colors <- mergedMEs$module_identity %>% unique()

mergedMEs$gene_ID_robin[mergedMEs$module_identity == colors[1]] 

# Make a list that returns gene names for a given cluster
module_color = colors
module = names(mergedMEs)
module_genes <- list()
module_color <- module_colors
# Get the genes from each of the modules
for (i in 1:length(module_color)) {
  module_genes[[i]] <- mergedMEs$gene_ID_robin[mergedMEs$module_identity == colors[i]] 
  names(module_genes)[[i]] <- module_color[[i]]
}



# load ejtk database

#install.packages("readxl")

library("readxl")


ejtk.db <- dbConnect(RSQLite::SQLite(), paste0(path_to_repo,"/data/databases/TC6_fungal_ejtk.db"))
DEGs <- read_excel('./supplement/upregulated_DEGs_Will_2020.xlsx') 


ll <- DEGs$arb2_gene


## MAKE YOUR LIST OF GENES OF INTEREST ##

# LIST ONE - WGCNA modules
list1 <- module_genes
sapply(list1, length)

## LIST TWO - rhythmic genes
list2 <- list(ll,ll)
names(list2) <- c('DEG_alive', 'copy')
sapply(list2, length)


## CHECK FOR OVERLAP

## make a GOM object
gom.1v2 <- newGOM(list1, list2,
                  genome.size = nGenes)
png(paste0(path_to_repo, "/results/figures/", sample.name,"_DEG_gom_1v2.png"), 
    width = 20, height = 30, units = "cm", res = 300)
drawHeatmap(gom.1v2,
            adj.p=T,
            cutoff=0.05,
            what="odds.ratio",
            # what="Jaccard",
            log.scale = T,
            note.col = "grey80")
trash <- dev.off()


############ DEG and rhy
# ejtk.db <- dbConnect(RSQLite::SQLite(), paste0(path_to_repo,"/data/databases/TC6_fungal_ejtk.db"))
dat2 <- read.csv(glue('{path_to_repo}/data/ophio_cflo_TC6_data.csv'))


# DEFINE GENES OF INTEREST

rhy.24 <- dat2 %>% 
  filter(GammaP_24h < 0.05) %>% 
  pull(gene_ID_robin) %>% 
  as.character()

rhy.12 <- dat2 %>% 
  filter(GammaP_12h < 0.05) %>% 
  pull(gene_ID_robin) %>% 
  as.character()

rhy.08 <- dat2 %>% 
  filter(GammaP_08h < 0.05) %>% 
  pull(gene_ID_robin) %>% 
  as.character()



## LIST TWO - rhythmic genes
list1 <- list(rhy.24, rhy.12, rhy.08)

names(list1) <- paste0(sample.name, c("-24h", "-12h", "-08h"))
sapply(list1, length)

gom.1v2 <- newGOM(list1, list2,
                  genome.size = nGenes)
png(paste0(path_to_repo, "/results/figures/", sample.name,"_DEG_rhy_gom_1v2.png"), 
    width = 20, height = 30, units = "cm", res = 300)
drawHeatmap(gom.1v2,
            adj.p=T,
            cutoff=0.05,
            what="odds.ratio",
            # what="Jaccard",
            log.scale = T,
            note.col = "grey80")
trash <- dev.off()


#### Overlap all rhythmic genes and DEG

# check lenght of all rhy genes = 2921 
# length(rhy.24) + length(rhy.12) + length(rhy.08)

list1 <- c(rhy.24, rhy.12, rhy.08) %>% 
  unique()

## LIST TWO - rhythmic genes
list1 <- list(list1, list1, list1)

names(list1) <- paste0(sample.name, c("all_rhy", "copy", 'copy'))
sapply(list1, length)

gom.1v2 <- newGOM(list1, list2,
                  genome.size = nGenes)
png(paste0(path_to_repo, "/results/figures/", sample.name,"_DEG_all_rhy_gom_1v2.png"), 
    width = 20, height = 30, units = "cm", res = 300)
drawHeatmap(gom.1v2,
            adj.p=T,
            cutoff=0.05,
            what="odds.ratio",
            # what="Jaccard",
            log.scale = T,
            note.col = "grey80")
trash <- dev.off()

