# Check SNPs included in EC-PRSs in csPRS
setwd('/path/to/csPRS/score_split_analysis/analysis_compare/ec_prs')
library(data.table)
library(stringr)

# 36 SNPs in EC-PRS
snp.ecprs <- fread('weight_ecPRS.txt',data.table=F)
### harmonize SNPs
#### UKB SNPs
snp.ukb <- fread('/path/to/data/genotypes/ukb_genotypes.bim',
		 data.table=F)
snp.ecprs.harm <- merge(snp.ecprs,snp.ukb,
			by.x=c('Chromosome','pos_hg19'),
			by.y=c('V1','V4'),
			all.x=T)
sum(snp.ecprs.harm$Variant != snp.ecprs.harm$V2) # every SNP is present!!!!!

# csPRS GWAS with cell type annotated
path='/path/to/csPRS/score_cs_beta/gwas_cs/'
files <- list.files(path,'.txt')
cell <- str_remove_all(files,'.txt')
snp.csprs <- c()
for(c in cell){
  tmp <- fread(paste0(path,c,'.txt'),data.table=F)
  tmp$cell <- c

  snp.csprs <- rbind(snp.csprs,tmp)
}

#### Add cell annotations for EC-PRS SNPs
setdiff(snp.ecprs$Variant,snp.csprs$SNP) # 26 missing......
snp.ecprs2 <- snp.ecprs[snp.ecprs$Variant %in% snp.csprs$SNP,] # 9 left
annot <- lapply(snp.ecprs2$Variant,function(x){
		  add <- snp.csprs$cell[snp.csprs$SNP==x]
		  out <- paste(unique(add),collapse=';')
			})
annot.res <- do.call(rbind,annot)
out <- data.frame(snp=snp.ecprs2$Variant,
		  cell=annot.res)

write.csv(out,'snp_ecprs_compare.csv',row.names=F)
