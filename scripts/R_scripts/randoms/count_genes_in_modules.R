set.seed(420)
rm(list = ls())

#' Load the libraries
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted, WGCNA, igraph)
pacman::p_load(patchwork)

#' set conflict preference
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
conflict_prefer("hclust", "flashClust")

#' set path to your working directory
path_to_repo = "."

# script name
script.name = ""

#' load functions
# customized theme for publication quality figures
source(paste0(path_to_repo,"/functions/theme_publication.R"))
# function to perform enrichment analysis
source(paste0(path_to_repo,"/functions/enrichment_analysis.R"))
# function to plot z-scores (Cflo genes only)
source(paste0(path_to_repo,"/functions/plot_zscores.R"))


pacman::p_load(GeneOverlap)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf
pacman::p_load(glue)


# Make a list that returns gene names for a given cluster
module_color = colors
module = names(mergedMEs)
module_colors <-
  data.frame(module_label=module) %>%
  mutate(module_color = str_replace(module_label, "ME", ""))


### Start getting overlap numbers -----------------------------------------------

species <- 'ophio_cflo'

### Data of genes with their modules
all.modules <- read_csv(glue('./results/networks/{species}_gene_IDs_and_module_identity.csv'))

# Get all the module colos
module_colors <- 
  all.modules %>% 
  select(module_identity) %>% 
  unique()

module_genes <- list()
module_color <- module_colors$module_color
# Get the genes from each of the modules
for (i in 1:length(module_colors)) {
  
  module_genes[[i]] <- names(all.modules)[which(module_colors==module_colors[[i]])]
  names(module_genes)[[i]] <- module_colors[[i]]
}
# check the result | works
names(module_genes)
module_genes['salmon']

writeLines("#####################################################
How many genes are in each of my geneset of interest?
#####################################################")

## MAKE YOUR LIST OF GENES OF INTEREST ##


# load ejtk database
ejtk.db <- dbConnect(RSQLite::SQLite(), paste0(path_to_repo,"/data/databases/TC6_fungal_ejtk.db"))

# DEFINE GENES OF INTEREST
sample.name <- 'ophio_cflo'

rhy.24 <-
  tbl(ejtk.db, paste0(sample.name,"_zscores_24h")) %>% 
  filter(GammaP < 0.05) %>% pull(ID)

rhy.12 <-
  tbl(ejtk.db, paste0(sample.name,"_zscores_12h")) %>% 
  filter(GammaP < 0.05) %>% pull(ID)

rhy.08 <- 
  tbl(ejtk.db, paste0(sample.name,"_zscores_08h")) %>% 
  filter(GammaP < 0.05) %>% pull(ID)


# DRGs


# LIST ONE - WGCNA modules
list1 <- module_genes
writeLines("List of interesting genes #1
----------------------------
Genes in each of the identified gene-clusters or modules")
sapply(list1, length)

## LIST TWO - rhythmic genes
list2 <- list(rhy.24, rhy.12, rhy.08)
writeLines("List of interesting genes #2
----------------------------
Rhythmic genes in control Cflo heads")
names(list2) <- paste0(sample.name, c("-24h", "-12h", "-08h"))
sapply(list2, length)


## CHECK FOR OVERLAP

## make a GOM object
gom.1v2 <- newGOM(list1, list2,
                  genome.size = nGenes)
png(paste0(path_to_repo, "/results/figures/", sample.name,"_gom_1v2.png"), 
    width = 20, height = 30, units = "cm", res = 300)
drawHeatmap(gom.1v2,
            adj.p=T,
            cutoff=0.05,
            what="odds.ratio",
            # what="Jaccard",
            log.scale = T,
            note.col = "grey80")
trash <- dev.off()

# writeLines("Visualizing the significant overlaps between your lists of interesting genes and the identified modules")


# 
# The above plot identifies the gene-clusters that contain most of the genes that show daily rhythms in healthy Cflo heads (Cflo-healthy-controls).
# 
# HOW TO READ THE FIGURE?
#   
#   - Briefly explain here
# - Talk about Odds-Ratio
# - Darker the green, more significant is the overlap (also indicated by the adj_pval)
# 
# Next, we can try to identify the ant gene-clusters that underlie behavioral plasticity, as well as the ant clusters that are affected during behavioral manipulation induced by *Ophiocordyceps*.
# 
# foo %>% 
#   mutate(GCN = "beau") %>% 
#   write.csv(paste0(path_to_repo,"/results/networks/", sample.name,"_genes_and_module_identity.csv"),
#             row.names = F)


### DO IT MY OWN ----------------------------------------------------------

species <- 'ophio_cflo'

### Data of genes with their modules
all.modules <- read_csv(glue('./results/networks/{species}_gene_IDs_and_module_identity.csv'))

# Get all the module colors in list
module_colors <- 
  all.modules %>% 
  pull(module_identity) %>% 
  unique()

## Give lenght of all modules
list_x <- NULL

for (i in module_colors) {
  module <- 
    all.modules %>% 
    filter(module_identity==i) %>% 
    pull(gene_ID_ncbi)
    
  len_module <- length(module)
  list_name <- glue('{i} ({len_module})')
  
  list_x <- append(list_x, list_name)
}

# write to txt
lapply(list_x, write, glue("./results/networks/{species}_length_module_list.txt"), append=TRUE, ncolumns=1000)

##### Get overlap numbers

module_name <- 'darkgrey'

module <- 
  all.modules %>% 
  filter(module_identity==module_name) %>% 
  pull(gene_ID_ncbi)

intersect(module, rhy.08) %>% length()

