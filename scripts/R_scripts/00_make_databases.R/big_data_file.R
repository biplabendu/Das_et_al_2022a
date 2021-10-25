### Create big data file #######################
# That contains all the data for the TC6 ophio and beau

# Housekeeping ---------------------------------------------------------------
#
set.seed(420)
rm(list = ls())
#
## Load packages --
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
pacman::p_load(glue)
#
# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
#
## set parameters and thresholds
#
# Specify your path
path <- "/Users/roos_brouns/Dropbox/Ant-fungus/02_scripts/Git_Das_folder2/Das_et_al_2022a"
#
# Specify species
species <- 'ophio_cflo'



# 1. FPKM table  ----------------------------------------------------------
#
# Load in the expression data
exp.val <- read.csv(glue("{path}/results/normalized_gene_exp/raw_fpkm/{species}/normalized_gene_exp_{species}_all_samples.csv"),
                    header = T, stringsAsFactors = F, na.strings = c(NA, "", " ")) %>%
  as.data.frame()

# also make a colum for the gene ID's from Robin
#
# Load in the gene IDs
gene.IDs <- read.csv(glue("{path}/data/input/{species}/{species}_gene_names_robin_ncbi.csv"),
                     header = T, stringsAsFactors = F, na.strings = c(NA, "", " ")) %>%
  as.data.frame() 

# Join the data
data <- left_join(exp.val, gene.IDs, by = 'gene_ID_ncbi')
# rearange cols
data <- data[,c(16,1:15)]

# A. Expressed genes ---------------------------------------------------------

# list of all beau "Expressed" genes (>1 FPKM during the 24h period)
expressed.list <-
  data %>% 
  select(-end,
         -start)%>% 
  collect() %>%
  na.omit() %>% # Remove NA's
  filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>% # expression > 1 FPKM 
  pull(gene_ID_ncbi) %>% # get the gene names of the expressed genes
  unique() # check wheter there are duplicates

# Make a df with a yes when the gene is expressed and no if not
expressed.df <-
  data %>%
  select(gene_ID_ncbi, gene_ID_robin) %>%
  mutate(expressed = ifelse(gene_ID_ncbi %in% expressed.list,"yes","no"))

# Join the expressed df to our data df
data <- left_join(data, expressed.df, by = c('gene_ID_ncbi', 'gene_ID_robin'))

# B. Rhythmic genes ---------------------------------------------------------
# Genes are rhythmic when GammaP > 0.05

periods <- c('08',12,24) %>%
  as.character

for (i in periods) {
  # name the file with the species and period
  name <- glue('{species}_zscores_{i}h')
  # get the CSV name
  csv.name <- glue("{path}/results/ejtk_output/{species}/zscore/{species}_zscores_noNAs_cos{i}_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt")
  #
  # Note: ophio.cflo.24.zscore is exactly the same as ophio.cflo.24 meaning that using
  #       FPKM values or z-scores to run eJTK does not change the result
  # 
  # load in the zscore data 
  zscore <- read.csv(csv.name,
                     sep = "\t", header = T, stringsAsFactors = F)
  
  # all genes with a GammaP < 0.05 are rhythmic and get a 'yes'
  rhythmic <-
    zscore %>%
    select(gene_ID_ncbi = ID, 
           GammaP) %>%
    collect %>%
    mutate(rhythmic = ifelse(GammaP < 0.05, "yes", "no"))

  # Rename the cols
  colnames(rhythmic)[2] <- glue('GammaP_{i}h')
  colnames(rhythmic)[3] <- glue('rhythmic_{i}h')
  
  # Join with data df
  data <- left_join(data, rhythmic, by = 'gene_ID_ncbi')

  }

# C. Identified modules ----------------------------------------------
#
# Load in data with assigned modules
 modules <- read.csv(glue('{path}/results/networks/{species}_gene_IDs_and_module_identity.csv')) %>%
  select(gene_ID_ncbi, module_identity)

data <- left_join(data, modules, by = 'gene_ID_ncbi')

# D. Rhythmic module -----------------------------------------
# If the module has significant overlap with rhytmic genes 

if (species == 'ophio_cflo') {
  # For ophio the modules with 24h_rhytmic overlap are tan, mighnightblue, brown, darkturquiose
  rhythmic_module <-
  data %>%
  select(gene_ID_ncbi, 
         module_identity) %>%
  collect %>%
  mutate(rhythmic_module = ifelse(module_identity == 'tan' | module_identity == 'mignightblue' | module_identity == 'brown' | module_identity == 'darkturquiose' , "yes", "no")) %>%
  select(-module_identity)
} else if (species == 'beau') {
  # For ophio the modules with 24h_rhytmic overlap are darkred, red, salmon
  rhythmic_module <-
    data %>%
    select(gene_ID_ncbi, 
           module_identity) %>%
    collect %>%
    mutate(rhythmic_module = ifelse(module_identity == 'red' | module_identity == 'darkred' | module_identity == 'salmon', "yes", "no")) %>%
    select(-module_identity)
} else {
  print('We do not have data for this species')
}

# Rename the cols
colnames(rhythmic_module)[2] <- glue('rhythmic_module_24h')

# Join with data df
data <- left_join(data, rhythmic_module, by = 'gene_ID_ncbi')

# E. Annotations --------------------------------------------

if (species == 'ophio_cflo') {
## Ophio_cflo
annots <- read.csv(glue('{path}/data/input/{species}/{species}_annots_robin_ncbi.csv')) %>%
  select(-c(arb2_gene, sc16a_homolog)) %>%
  mutate(signalP = ifelse(!is.na(annots$signalP), 'yes', 'no')) %>%
  mutate(SSP = ifelse(!is.na(annots$SSP), 'yes', 'no')) %>%
  mutate(TMHMM = ifelse (!is.na(annots$TMHMM), 'yes', 'no'))
} else if (species == 'beau') {
## Beau
annots <- read.csv(glue('{path}/data/input/{species}/{species}_annots_robin_ncbi.csv')) %>%
  select(-c(gene_ID_robin, Start, End)) %>%
  mutate_all(list(~na_if(.,""))) %>%
  mutate(signalP = ifelse(!is.na(annots$signalP), 'yes', 'no')) %>%
  mutate(SSP = ifelse(!is.na(annots$SSP), 'yes', 'no')) %>%
  mutate(TMHMM = ifelse (!is.na(annots$TMHMM), 'yes', 'no'))
} else {
  print("We do not have data for this species")
}

annots <- read.csv(glue('{path}/data/input/{species}/{species}_annots_robin_ncbi.csv'))  %>%
  select(-c(gene_ID_robin, Start, End)) %>%
  mutate_all(list(~na_if(.,""))) %>%
  mutate(signalP = ifelse(is.na(annots$signalP), 'no', 'yes')) %>%
  mutate(SSP = ifelse(is.na(annots$SSP), 'no', 'yes')) %>%
  mutate(TMHMM = ifelse (!is.na(annots$TMHMM), 'yes', 'no'))


data <- left_join(data, annots, by = 'gene_ID_ncbi')

# Othologs ----------------------------------------------

if (species == 'ophio_cflo'){
  ## Ophio_cflo
  orthos <- read.csv(glue('{path}/data/input/{species}/{species}_annots_robin_ncbi.csv')) %>%
    select(c(gene_ID_ncbi, ophio_kim_homolog = sc16a_homolog))
} if else (species == 'beau') {
  ## Beau
  print ('No ortholog available for Beau')
} else {
  print("We do not have data for this species")
}

data <- left_join(data, orthos, by = 'gene_ID_ncbi')

#### Write to csv file
write.csv(data, glue("{path}/data/{species}_TC6_data.csv"), row.names = FALSE)





