
# 01. Housekeeping ----------------------------------------------------------
#
## Load libraries
library(tidyverse)
#
## Set parameters and thresholding
gamma.pval = 0.05
#
## End.

# 02. Ophiocordyceps camponoti-floridani --------------------------------------
#
## Period = 24h
ophio.cflo.24 <- read.csv("./results/ejtk_output/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
           sep = "\t", header = T, stringsAsFactors = F)
ophio.cflo.24 <- ophio.cflo.24 %>%  
  filter(GammaP < gamma.pval) %>% 
  pull(ID)
#
#
## Period = 12h
ophio.cflo.12 <- read.csv("./results/ejtk_output/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples_cos12_ph0020by4_as0420by4_jtkout_GammaP.txt",
                          sep = "\t", header = T, stringsAsFactors = F)
ophio.cflo.12 <- ophio.cflo.12 %>% 
  filter(GammaP < gamma.pval) %>% 
  pull(ID)
#
#
## Period = 8h
ophio.cflo.08 <- read.csv("./results/ejtk_output/ophio_cflo/normalized_gene_exp_ophio_cflo_all_samples_cos8_ph0020by4_as0420by4_jtkout_GammaP.txt",
                          sep = "\t", header = T, stringsAsFactors = F)
ophio.cflo.08 <- ophio.cflo.08 %>% 
  filter(GammaP < gamma.pval) %>% 
  pull(ID)
#
#
## Save the rhythmic gene-lists as a data file
save(ophio.cflo.08,
     ophio.cflo.12,
     ophio.cflo.24,
     file = "./data/rdata/rhythmic_genes/ophio_cflo_rhy_genes_gammap_05.RData")
#
## To read the RData file, use:
## load("/path/to/file.RData)
## End.

# 03. Ophiocordyceps kimflemingae ---------------------------------------------
#
## Light-Dark cycle
  ## Period = 24h
ophio.kim.ld.24 <- read.csv("./results/ejtk_output/ophio_kim/ld/normalized_gene_exp_ophio_kim_ld_samples_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
ophio.kim.ld.24 <- ophio.kim.ld.24 %>% 
  filter(GammaP < gamma.pval) %>% 
  pull(ID)
#
#
## Dark-Dark
  ## Period = 24h
ophio.kim.dd.24 <- read.csv("./results/ejtk_output/ophio_kim/dd/normalized_gene_exp_ophio_kim_dd_samples_cos24_ph0020by4_as0420by4_jtkout_GammaP.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
ophio.kim.dd.24 <- ophio.kim.dd.24 %>% 
  filter(GammaP < gamma.pval) %>% 
  pull(ID)
#
#
## Save the rhythmic gene-lists as a data file
save(ophio.kim.dd.24,
     ophio.kim.ld.24,
     file = "./data/rdata/rhythmic_genes/ophio_kim_rhy_genes_gammap_05.RData")
#
##
## End.