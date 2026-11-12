# generate jobs for SDPR
setwd('/path/to/csPRS/effect_size_comparison/code')
library(data.table)
library(stringr)

for(c in c('type','state')){
  path <- paste0('/path/to/csPRS/score/weights/weights_specificity_top10pct_cell_',c)
  files <- list.files(path,'txt')
  prefix <- str_remove_all(files,'.txt')
  prefix <- setdiff(prefix,c("Atrial-Cardiomyocyte.txt","Ventricular-Cardiomyocyte.txt"))
  if(c=='state'){
    prefix <- setdiff(prefix,c('Meso','doublets','nan')) # remove duplicates
  }
  for(p in prefix){
    #### Prepare GWAS to have only cs variants
    gwas <- fread(paste0(path,'/',p,'.txt'),data.table=F)
    gwas_clean <- subset(gwas,select=c('SNP','EA','NEA','BETA','P')) 
    colnames(gwas_clean) <- c('SNP','A1','A2','BETA','P')

    write.table(gwas_clean,paste0('../gwas_clean/gwas_',p,'.txt'),
		row.names=F,sep='\t',quote=F)
  }
}

# Generate jobs for SDPR
samplesize=184305
path_to_gwas <- '/path/to/csPRS/effect_size_comparison/gwas_clean/'
gwas <- list.files(path_to_gwas,'.txt')
gwas <- setdiff(gwas,c('gwas_Atrial-Cardiomyocyte.txt','gwas_Ventricular-Cardiomyocyte.txt'))
prefix <- str_remove_all(gwas,'gwas_')
prefix <- str_remove_all(prefix,'.txt')
path_out='/path/to/csPRS/effect_size_comparison/cs_weight/by_chr/'
job <- c()
for(p in prefix){
  for(chr in 1:22){
    add <- paste0('cd /path/to/software/SDPR/;./SDPR -mcmc -ref_dir /path/to/resources/sdpr_reference/ -ss ',
  		  path_to_gwas,'gwas_',p,'.txt -N ',samplesize,' -chr ',chr,
		  ' -out ',path_out,p,'_chr',chr,' -n_threads 3')
    job <- c(job,add)
   }
}

write.table(job,'04_cell_sdpr.job',row.name=F,col.names=F,sep='',quote=F)

