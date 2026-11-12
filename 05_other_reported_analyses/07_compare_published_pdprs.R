# Compare the csPRSs to PD-PRSs decomposed before
setwd('/path/to/csPRS/score_split_analysis/analysis_compare/pd_prs')
library(data.table)

# Read csPRSs
path.csprs <- '/path/to/csPRS/score_split_analysis/data/'
dat.csprs <- fread(paste0(path.csprs,'ukb_clustercsprs.txt.gz'),data.table=F)
### extract csPRSs
csprs <- colnames(dat.csprs)[grep('csbeta',colnames(dat.csprs))]
dat.csprs2 <- subset(dat.csprs,select=c('eid',csprs))

# CAD phenotype and batches
dat.cad <- fread(paste0(path.csprs,'ukb_inciCAD.txt.gz'),data.table=F)

dat <- merge(dat.cad,dat.csprs2,by='eid')
### remove training
dat <- dat[dat$batch!='training',]

# Read pdPRSs
path.pdprs <- '/path/to/resources/published_scores/'
files <- list.files(path.pdprs,'profile')
prefix <- c('basic_condition','bp','cvd','immune_system','lipids',
	    'others','obesity','respiratory','t2d')
### read one by one and merge
for(i in 1:length(files)){
  tmp <- fread(paste0(path.pdprs,files[i]),select=c('IID','SCORE'))
  colnames(tmp) <- c('eid',paste0('pdprs_',prefix[i]))

  if(i==1){
    dat.pdprs <- tmp
  }else{
    dat.pdprs <- merge(dat.pdprs,tmp,by='eid')
  }
}
pdprs <- colnames(dat.pdprs)[-1]

# Merge with csPRSs and CAD
dat <- merge(dat,dat.pdprs,by='eid')
dim(dat)

### Correlations between csPRSs and pdPRSs
library(Hmisc)
library(tidyverse)

res <- rcorr(
  as.matrix(dat[, c(csprs, pdprs)]),
  type = "pearson"
)

# tidy output (csPRS vs PD-PRS only)
cor_df <- expand_grid(
  var1 = csprs,
  var2 = pdprs
) %>%
  mutate(
    cor = map2_dbl(var1, var2, ~ res$r[.x, .y]),
    p   = map2_dbl(var1, var2, ~ res$P[.x, .y])
  )

write.csv(cor_df,'corr_pdprs.csv',row.names=F)
