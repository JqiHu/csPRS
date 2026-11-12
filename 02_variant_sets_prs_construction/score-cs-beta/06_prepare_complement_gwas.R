# Prepare remaining cs-gwas subsets for SDPR as input
setwd('/path/to/csPRS/score_cs_beta/gwas_cs_remain')
library(data.table)
library(stringr)

# Read SDPR BETA GWAS
gwas <- fread('/path/to/resources/prs/sumstats/clean_for_SDPR/SDPR_cad_cardiogramc4a_sum_rm_dup_overlapUKB.txt',
	      data.table=F)

## csbeta
path='../gwas_cs/'
files <- list.files(path,'.txt')
cells <- str_remove_all(files,'.txt')

# Run group by group
num.snp <- c()
for(c in cells){
  csgwas <- fread(paste0(path,c,'.txt'),data.table=F)
  out <- gwas[which(gwas$SNP %in% csgwas$SNP==F),]
  print(nrow(out))

  out <- out[duplicated(out$SNP)==F,]
  write.table(out,paste0(c,'.txt'),
	      row.names=F,sep='\t',quote=F)

  num.snp <- c(num.snp,nrow(out))
}

num.snp.out <- data.frame(group=cells,
			  num=num.snp)
write.csv(num.snp.out,'num_snp.csv',row.names=F)


