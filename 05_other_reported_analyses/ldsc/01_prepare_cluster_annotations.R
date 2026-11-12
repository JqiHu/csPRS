# Generate new annotations based on gene expression in cells
setwd('/path/to/csPRS/heritability_par/annot_specificity_top10pct_merge251127')
library(data.table)
library(stringr)

### SNPs assigned to Genes expressed in specific cells
cell_type <- fread('/path/to/csPRS/score/weights/weights_by_cell_type/SDPR_CAD_cell_type.txt.gz',data.table=F)

cell_state <- readRDS('/path/to/csPRS/score/weights/weights_by_cell_state/SDPR_CAD_cell_state.rds')
cell_state2 <- do.call(rbind,cell_state)

snp.all <- rbind(cell_type,cell_state2)
snp.all$cell <- str_replace_all(snp.all$cell,'-','_')

### 1000G from LDSC
path='/path/to/software/ldsc/ref/1000G_EUR_Phase3_plink/'
ref <- c()
for(chr in 1:22){
  add <- fread(paste0(path,'1000G.EUR.QC.',chr,'.bim'))
  ref <- rbind(ref,add)
}
ref$index <- paste0(ref$V1,':',ref$V4)

#### Function to write out annotations chr by chr
getOut <- function(file,prefix){
  for(chr in 1:22){
    out <- file[file$CHR==chr,]
    write.table(out,gzfile(paste0('Annot_',prefix,'_chr',chr,'.annot.gz')),
                row.names=F,sep=' ',quote=F)
    print(paste0('Finished for chr',chr))
  }
}

### Create annotations across cell clusters
###### Read genes
info <- readxl::read_xlsx('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127/jaccard_similarity_group_edit.xlsx')
genes <- readRDS('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127/top10_group_10.rds')
## check
(setdiff(unique(snp.all$cell),unique(info$Var1)))

res <- data.frame(CHR=ref$V1,
		  BP=ref$V4,
		  SNP=ref$V2,
		  CM=ref$V3)

for(c in unique(info$Cluster2_name)){
  cells <- info$Var1[info$Cluster2_name==c]  
  tmp <- snp.all[snp.all$cell %in% cells,]
  tmp <- tmp[duplicated(tmp$SNP)==F,]
  genes.tmp <- genes[[c]]
  ### extract SNPs
  tmp.top <- tmp[tmp$gene_name %in% genes.tmp,]
  snps <- ref$index %in% tmp.top$index
  if(sum(snps)==0){
    print(paste0('Skip ', c))
    next
  }

  annot <- ifelse(snps,1,0)
  print(table(annot))
  res <- cbind.data.frame(res,annot)
  colnames(res)[ncol(res)] <- c
}
getOut(res,'merge251127')
