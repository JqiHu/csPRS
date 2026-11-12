# Prepare cs-gwas subsets for SDPR as input
setwd('/path/to/csPRS/score_cs_beta/gwas_cs')
library(data.table)
library(stringr)

# Read cell cluster info
info <- readxl::read_xlsx('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127/jaccard_similarity_group_edit.xlsx')

# cleaned GWAS for top 10% specificity 
path='/path/to/csPRS/effect_size_comparison/gwas_clean'
files <- list.files(path)
cells <- str_remove_all(files,'gwas_')
cells <- str_remove_all(cells,'.txt')

# check
(setdiff(cells,unique(info$Var1)))
info$Var1 <- str_replace_all(info$Var1,'_','-')
(setdiff(cells,unique(info$Var1)))

# Run group by group
num.snp <- c()
for(g in unique(info$Cluster2_name)){
  tmp.cells <- info$Var1[info$Cluster2_name==g]
  out <- c()
  for(c in tmp.cells){
    gwas <- fread(paste0(path,'/gwas_',c,'.txt'),data.table=F)  
    out <- rbind(out,gwas)
  }
  # remove duplicates
  out <- out[!duplicated(out$SNP), , drop = FALSE]

  write.table(out,paste0(g,'.txt'),
	      row.names=F,sep='\t',quote=F)

  num.snp <- c(num.snp,nrow(out))
}

num.snp.out <- data.frame(group=unique(info$Cluster2_name),
			  num=num.snp)

write.csv(num.snp.out,'num_snp.csv',row.names=F)
