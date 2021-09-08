
### Add the ncbi gene names to the ophio_cflo annotation file
#
foo <- read.csv("./data/input/ophio_cflo/FullBlast_EC05_RNAseq_orignal_copy_26Aug19.csv",
                header = T,
                stringsAsFactors = F)

bar <- read.csv("./data/input/ophio_cflo/gffs/ophio_kim_gene_names_robin_ncbi.csv",
                header = T,
                stringsAsFactors = F)

bar %>% head()

foo %>% 
  left_join(bar[,3:4], by=c("arb2_gene" = "attributes_robin")) %>% 
  select(arb2_gene,
         gene_name_ncbi = attributes_ncbi,
         blast_annot = Description.y,
         GOs = GO,
         pfams = PFAM,
         signalP = signalp,
         SSP,
         TMHMM,
         sc16a_homolog=sc16a_gene) %>%
  unique() %>% 
  write.csv(.,
            file = "./data/input/ophio_cflo/ophio_cflo_annots_robin_ncbi.csv",
            row.names = F)
