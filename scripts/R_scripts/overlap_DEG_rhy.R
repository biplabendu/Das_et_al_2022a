### What is in the annoataions of overlapping genes
set.seed(420)
rm(list = ls())

pacman::p_load(pheatmap, dendextend, tidyverse, viridis, ggthemes)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
pacman::p_load(glue)


# Read in data
path = "/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a"
overlap_data <- read_csv(glue('{path}/supplement/overlapping_genes_with_DEG.csv'))
annot_data <- read_csv(glue('{path}/data/input/ophio_cflo/complete_annotations/FullBlast_EC05_RNAseq_orignal_copy_26Aug19.csv')) %>% 
  select()

# 24h overlap with upregulated DEG TMHMM
DEG.TMHMM <- 
  overlap_data$`124 common elements in DEG TMHMM and 24h rhythmic` 

DEG.signalP <- 
  overlap_data$`55 common elements in DEG signalP and 24h rhytmic`

DEG.module_A7 <- 
  overlap_data$`127 common elements in DEGS and Modele A7`

DEG.module_A1 <- 
  overlap_data$`174 common elements in DEGS and Modele A1`

overlap.w.annots <- 
  annot_data %>% 
  mutate(rhy24h_DEG_TMHMM_overlap = ifelse(arb2_gene %in% DEG.TMHMM, 'yes', 'no')) %>% 
  mutate(rhy24h_DEG_signalP_overlap = ifelse(arb2_gene %in% DEG.signalP, 'yes', 'no')) %>% 
  mutate(rhy24h_DEG_mA7_overlap = ifelse(arb2_gene %in% DEG.module_A7, 'yes', 'no')) %>% 
  mutate(rhy24h_DEG_mA1_overlap = ifelse(arb2_gene %in% DEG.module_A1, 'yes', 'no'))
## Need one last filter step to filter the 4 cols with all no out. 


write_csv(overlap.w.annots, glue('{path}/supplement/overlap_24hrhy_upDEG_w_annots_filtered.csv'))
