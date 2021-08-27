
### Cleaning the Ophio_kim annotation file
#
# foo <- read.csv("./results/normalized_gene_exp/ophio_kim/ophio_kim_annots_robin.csv",
#                 header = T,
#                 stringsAsFactors = F)
# 
# library(tidyverse)
# 
# protein_name <- character()
# for (i in 1:nrow(foo)) {
# protein_name[i] <- str_split(foo$ProteinId.t1, ".t1")[[i]][1]
# }
# 
# foo[1] <- protein_name
# 
# foo %>% 
#   select(gene_name=ProteinId.t1,
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
#
#### Done.

### Add the ncbi gene names to the annotation file
# foo <- read.csv("./results/normalized_gene_exp/ophio_kim/ophio_kim_annots_robin_cleaned.csv",
#                 header = T, stringsAsFactors = F)
# bar <- read.csv("./results/normalized_gene_exp/ophio_kim/ophio_kim_gene_names_robin_ncbi.csv",
#                 header = T, stringsAsFactors = F)
# 
# foo %>% 
#   left_join(bar[,(3:4)], by=c("gene_name" = "attributes_robin")) %>% 
#   select(gene_name,
#          gene_name_ncbi=attributes_ncbi,
#          everything()) %>% 
#   write.csv(.,
#             file = "./results/normalized_gene_exp/ophio_kim/ophio_kim_annots_robin_ncbi.csv",
#             row.names = F)
#
#### Done.
