# Merge certain cells based on Jaccard similarity
setwd('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127')
library(readxl)
library(data.table)

# read groups defined
group <- read_xlsx('jaccard_similarity_group_edit.xlsx')

## Read all 77 gene lists
gene1 <- readRDS('../cell_type/top10_gene.rds')
gene1 <- gene1[-which(names(gene1)=='doublets')]
gene2 <- readRDS('../cell_state/top10_gene.rds')

genes <- c(gene1,gene2)

# Create new list for each group
out <- list()
for(g in unique(group$Cluster2_name)){
  tmp.cell <- group$Var1[group$Cluster2_name==g]
  tmp.gene <- unlist(genes[tmp.cell])
  tmp.gene2 <- unique(tmp.gene)

  out[[g]] <- tmp.gene2
}
saveRDS(out,'top10_group_10.rds')

### number of genes 
num.out <- data.frame(group=names(out),
		      num=lengths(out))

write.csv(num.out,'top10_group10_num_gene.csv',row.names=F)
