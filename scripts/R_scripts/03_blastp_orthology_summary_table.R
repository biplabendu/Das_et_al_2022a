# setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/feature_tables")

setwd("~/Documents/GitHub/Das_et_al_2022a/data/proteinortho/")

rm(list = ls())

# Load libraries
pacman::p_load(conflicted, tidyverse)

# Set preference for functions
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("setdiff", "BiocGenerics")
conflict_prefer("union", "BiocGenerics")
conflict_prefer("intersect", "BiocGenerics")

source("./../../functions/theme_publication.R")

# Read the results of proteinortho ----------------------------------------

dat <- read.csv("./protein_fasta/02_results/ophcf2_v_beau_1.csv",
                 header = T, stringsAsFactors = F, na.strings = c(""," ","NA"))


# Clean the data ----------------------------------------------------------

## Keep only the one-to-one orthologs

# which ophio proteins have multiple beau orthologs?
duplicate.ophio.proteins <- 
  dat %>% 
  group_by(ophcf2_protein) %>% 
  summarize(num = n()) %>% 
  filter(num>1) %>% 
  pull(ophcf2_protein)

# which beau proteins have multiple ophio orthologs?
duplicate.beau.proteins <- 
  dat %>% 
  group_by(beau_protein) %>% 
  summarise(num = n()) %>% 
  filter(num>1) %>% 
  pull(beau_protein)

# remove the ophio and beau proteins that have multiple ortholog hits
dat.filtered <- 
  dat %>% 
  # remove the ophio proteins that have multiple hits for beau orthologs, and vice versa
  filter(!ophcf2_protein %in% duplicate.ophio.proteins) %>% 
  filter(!beau_protein %in% duplicate.beau.proteins)
  
# # double check to see if it worked (it did!)
# dat.filtered %>% 
#   group_by(ophcf2_protein) %>% 
#   # group_by(beau_protein) %>% 
#   summarize(num = n()) %>% 
#   filter(num>1)
  
# check to see the distribution of bitscore (or evalue) makes sense
dat.filtered %>% 
  ggplot(aes(log2(bitscore_ab))) +
  geom_vline(xintercept = log2(50), col="red") + # usually a bitscore ≥50 indicates sig. homology: https://www.biostars.org/p/187230/
  geom_histogram() +
  theme_Publication()
  

# Get the gene_names for the protein IDs ----------------------------------

## The following two csvs were downloaded from NCBI Genome website
    # Download genome annotation > Tabular format

## Read the files that have protein name to gene name mappings
# OPHIO
ophio.meta <- 
  read.csv("./feature_tables/ophcf2_protein_gene_names_ncbi.csv", 
                       header = T, stringsAsFactors = F, na.strings = c(""," ","NA")) %>% 
  select(ophcf2_protein = Protein.product, ophio_gene = Locus.tag, ophio_desc = Protein.Name) %>% 
  distinct()
# BEAU
beau.meta <- 
  read.csv("./feature_tables/beau_protein_gene_names_ncbi.csv",
                      header = T, stringsAsFactors = F, na.strings = c(""," ","NA")) %>% 
  select(beau_protein = Protein.product, beau_gene = Locus.tag, beau_desc = Protein.Name) %>% 
  distinct()


## Combine gene_name info from the above two files into the proteinortho results (dat.filtered)
ophio.beau.homology <- 
  dat.filtered %>% 
  left_join(ophio.meta, by="ophcf2_protein") %>% 
  left_join(beau.meta, by="beau_protein") %>% 
  select(ophio_gene, beau_gene, ophio_desc, beau_desc, everything())

## Save the file for further reference
ophio.beau.homology %>% 
  write.csv("./../../results/proteinortho/ophio_beau_homology.csv", row.names = F)

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
