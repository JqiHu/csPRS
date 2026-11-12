# Analyze genes shared across cells
setwd('/path/to/csPRS/gene_list')
library(data.table)
library(stringr)
library(vegan)
library(proxy)
library(tidyverse)

# Type of genes
prefix='specitificy_top10'
#prefix='top500'
#prefix='top500qc'

## read all list
type <- readRDS(paste0('genes_',prefix,'/cell_type/top10_gene.rds'))
state <- readRDS(paste0('genes_',prefix,'/cell_state/top10_gene.rds'))
genes <- c(type,state)
#genes <- fread('genes_top500_select/genes_top500_all.csv',data.table=F)
### remove duplicated ones
#genes2 <- genes[-c(3,56,60)]
#files <- list.files('genes_top500_qc')
#genes <- c()
#for(i in 1:length(files)){
#  tmp <- fread(paste0('genes_top500_qc/',files[i]),data.table=F)
#  cell <- str_remove_all(files[i],'_top500qc.txt')
#  tmp$cell <- cell
#  genes <- rbind(genes,tmp)  
#}
genes$index <- paste0(genes$hgnc_symbol,':',genes$cell)
genes2 <- genes[duplicated(genes$index)==F,]
#genes2 <- genes2[-which(genes2$cell %in% c('Meso','nan')),]
genes2 <- split(genes2$hgnc_symbol, genes2$cell)

##### calculate the shared proportion
# Pairwise jaccard similarity
n <- length(genes2)
J <- matrix(NA_real_, n, n, dimnames = list(names(genes2), names(genes2)))

for (i in seq_len(n)) {
  for (j in i:n) {
    inter <- length(intersect(genes2[[i]], genes2[[j]]))
    uni   <- length(union(genes2[[i]], genes2[[j]]))
    J[i, j] <- J[j, i] <- if (uni == 0) 0 else inter / uni  # define 0 if both empty
  }
}

sim_jaccard <- J

# Inspect result
dim(sim_jaccard)
sim_jaccard[1:5, 1:5]

# Convert it to long matrix
jac_long <- as.data.frame(sim_jaccard) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Jaccard")

write.csv(jac_long,paste0('./shared_gene/',prefix,'.csv'),row.names=F)
