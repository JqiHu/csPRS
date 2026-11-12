# Organize PRSs
setwd('/path/to/csPRS/score_split_analysis/data')
library(data.table)
library(stringr)

# Path to score files
path='/path/to/csPRS/score_cs_beta/score/'
### files
files <- list.files(path,'.sscore')
cluster <- str_remove_all(files,'ukb_')
cluster <- str_remove_all(cluster,'.sscore')
#### read scores
for(i in 1:length(files)){
  tmp <- fread(paste0(path,files[i]),data.table=F,
	       select=c('IID','SCORE1_AVG','SCORE2_AVG'))
  colnames(tmp) <- c('eid',paste0(cluster[i],'_csbeta'),
		     paste0(cluster[i],'_beta'))
  if(i==1){
    dat <- tmp
  }else{
    dat <- merge(dat,tmp,by='eid')
  }
}
print(dim(dat))

write.table(dat,gzfile('ukb_clustercsprs.txt.gz'),
	    row.names=F,sep='\t',quote=F)

# Path to cs-remaining PRSs
path='/path/to/csPRS/score_cs_beta/score_remain/'
### files
files <- list.files(path,'.sscore')
cluster <- str_remove_all(files,'ukb_')
cluster <- str_remove_all(cluster,'.sscore')
#### read scores
for(i in 1:length(files)){
  tmp <- fread(paste0(path,files[i]),data.table=F,
               select=c('IID','SCORE1_AVG','SCORE2_AVG'))
  colnames(tmp) <- c('eid',paste0(cluster[i],'_csbeta'),
                     paste0(cluster[i],'_beta'))
  if(i==1){
    dat <- tmp
  }else{
    dat <- merge(dat,tmp,by='eid')
  }
}
print(dim(dat))

write.table(dat,gzfile('ukb_clustercsprs_remain.txt.gz'),
            row.names=F,sep='\t',quote=F)

