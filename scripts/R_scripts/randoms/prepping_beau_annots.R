## Cleaning the Beau annotation file >> NOT DONE

# foo <- read.csv('./Das_et_al_2022a/data/input/beau/Beaba1.summary.txt',
#                 header = T,
#                 sep = '\t',
#                 stringsAsFactors = F)

library(tidyverse)

# protein_name <- character()
# for (i in 1:nrow(foo)) {
# protein_name[i] <- str_split(foo$ProteinID.t1, ".t1")[[i]][1]
# }
# 
# foo[1] <- protein_name
# 
# foo %>%
#   select(gene_name=ProteinID.t1,
#          GOs=GO.annotation..terminal.nodes.only.,
#          pfams=PFAM.annotation,
#          signalP=SignalP.SmallSecretedProteins,
#          TMHMM=Transmembrane.domains,
#          sc16a_ccl_up_down=SC16a....CcL..Up.Down,
#          sc16a_ccd_up_down=SC16a....CcD..Up.Down,
#          ccl_ccd_up_down=CcL....CcD..Up.Down) %>%
#   write.csv(.,
#             file = "./results/normalized_gene_exp/ophio_kim/ophio_kim_annots_robin_cleaned.csv",
#             row.names = F)

### Done.

## Add the ncbi gene names to the annotation file
foo <- read.csv("./Das_et_al_2022a/data/input/beau/Beaba1.summary.txt",
                sep = '\t', header = T, stringsAsFactors = F)

bar <- read.csv("./Das_et_al_2022a/data/input/beau/beau_gene_names_robin_ncbi.csv",
                header = T, stringsAsFactors = F)

foo %>% 
  left_join(bar[,(3:4)], by=c("ProteinID" = "attributes_robin")) %>% 
  select(ProteinID,
         gene_name_ncbi=attributes_ncbi,
         everything()) %>% 
  write.csv(.,
            file = "./Das_et_al_2022a/functions/funcdata/beau_annots_robin_ncbi.csv",
            row.names = F)
