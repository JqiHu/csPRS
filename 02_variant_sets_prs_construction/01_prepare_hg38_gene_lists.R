# Generate genes expressed in specific cells
# and map to positions (hg38)
# !!!! R 4.4
setwd('/path/to/csPRS/gene_list')
library(data.table)
library(biomaRt) # for gene position

# First, read average expression per cell
dat1 <- fread('average_expression_by_celltype.csv',data.table=F)
dat2 <- fread('average_expression_by_cellstate.csv',data.table=F)
## Get gene list
genes <- union(dat1$gene,dat2$gene)
length(genes) # 33,538 

# ---------------
# Get gene ranges
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"   # try "uswest" or "asia" if needed
)
ranges <- getBM(filters="hgnc_symbol",
                attributes=c('ensembl_gene_id','hgnc_symbol',"entrezgene_id",
                             'chromosome_name','start_position','end_position'),
                values=genes,
                mart=ensembl)
print(dim(ranges)) # position in hg38; 26,513 genes
### clean
ranges <- ranges[ranges$chromosome_name %in% 1:22,]
ranges <- ranges[duplicated(ranges$ensembl_gene_id)==F,]
print(dim(ranges)) # position in hg38 # 21,402
### check missing genes
missing <- setdiff(genes,ranges$hgnc_symbol)
length(missing) # 12,148 genes missing

## output reference table
write.csv(ranges,'gene_info.csv',row.names=F)
write.table(missing,'gene_missing.txt',
	    row.names=F,sep='\t',quote=F)


# ---------------
# Get cell-specific genes [genes which non-zero expression]
#### cell types
cells <- colnames(dat1)[-1]
for(c in cells){
  tmp <- subset(dat1,select=c('gene',c))
  colnames(tmp)[2] <- 'type'

  out <- tmp$gene[tmp$type!=0]
  print(length(out))

  res <- ranges[ranges$hgnc_symbol %in% out,]
  print(dim(res))

  write.table(res,paste0('./genes_by_cell_type/',c,'_genes.txt'),
	      row.names=F,sep='\t',quote=F)
}
#### cell states
cells <- colnames(dat2)[-1]
for(c in cells){
  tmp <- subset(dat2,select=c('gene',c))
  colnames(tmp)[2] <- 'type'

  out <- tmp$gene[tmp$type!=0]
  print(length(out))

  res <- ranges[ranges$hgnc_symbol %in% out,]
  print(dim(res))

  write.table(res,paste0('./genes_by_cell_state/',c,'_genes.txt'),
              row.names=F,sep='\t',quote=F)
}
