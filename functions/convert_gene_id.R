#### Convert Gene ID ----------------------------------
# his function is made to convert the gene names of NCBI to the gene names of Robin.

covert_id <- function(kind, id.list) {
  if (kind == ncbi) {
    ## code to covert ncbi to robin names
    return(robin.id.list)
  }
  else if (kind == robin) {
    ## code to convert robin to ncbi
    return(ncbi.id.list)
  }
  else {
    print("This function can only convert NCBI gene ID's to Robin's gene ID's and vice versa. Please specify which kind you want to concert 'ncbi' or 'robin'. ")
  }
  
}

covert_id_ncbi_to_robin_beau <- function(ncbi.id.list) {
  ## code to convert ncbi gene ID's to Robins gene ID's 
  robin.id.list <- ## something
  return(robin.id.list)
}

# Read in the data
matching.ids.beau <- read.csv("./data/input/beau/beau_gene_names_robin_ncbi.csv", sep = ',', header = TRUE, row.names = 4)
ncbi.id.df <- read.delim('./results/beau_not_expressed_list.txt', sep = ' ', header = FALSE)

### Made file into list
#
# Read in the data
ncbi.id.file <- scan("./results/beau_not_expressed_list.txt", what="", sep="\n")
# Separate elements by one or more whitepace
ncbi.id.list <- strsplit(ncbi.id.file, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(ncbi.id.list) <- sapply(ncbi.id.list, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
ncbi.id.list <- lapply(ncbi.id.list, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above


# Find matching rownames from list/df and extract the attributes_robin to obtain Robin's gene ID's
q <- matching.ids.beau[rownames(matching.ids.beau) %in% ncbi.id.df$V1, ]
p <- list(q$attributes_robin)

r <- matching.ids.beau[rownames(matching.ids.beau) %in% ncbi.id.list, ]

