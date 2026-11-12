# Analyze genes shared across cells
setwd('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127')
library(data.table)
library(stringr)
library(vegan)
library(proxy)
library(tidyverse)

## read all list
genes <- readRDS('top10_group_10.rds')

##### calculate the shared proportion
# Pairwise jaccard similarity
n <- length(genes)
J <- matrix(NA_real_, n, n, dimnames = list(names(genes), names(genes)))

for (i in seq_len(n)) {
  for (j in i:n) {
    inter <- length(intersect(genes[[i]], genes[[j]]))
    uni   <- length(union(genes[[i]], genes[[j]]))
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

write.csv(jac_long,'jaccard_sim.csv',row.names=F)
