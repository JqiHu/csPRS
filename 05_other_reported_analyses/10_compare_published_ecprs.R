# Compare the csPRSs to EC PRS
setwd('/path/to/csPRS/score_split_analysis/analysis_compare/ec_prs')
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

# Read EC-PRS
ecprs <- fread('ukb_ecprs.profile',data.table=F)
ecprs2 <- ecprs[,c('IID','SCORE')]
colnames(ecprs2) <- c('eid','ecprs')

# Merge with csPRSs and CAD
dat <- merge(dat,ecprs2,by='eid')
dim(dat)

### Correlations between csPRSs and pdPRSs
library(Hmisc)
library(tidyverse)

res <- rcorr(
  as.matrix(dat[, c(csprs, 'ecprs')]),
  type = "pearson"
)

# tidy output (csPRS vs PD-PRS only)
cor_df <- expand_grid(
  var1 = csprs,
  var2 = 'ecprs'
) %>%
  mutate(
    cor = map2_dbl(var1, var2, ~ res$r[.x, .y]),
    p   = map2_dbl(var1, var2, ~ res$P[.x, .y])
  )

write.csv(cor_df,'corr_ecprs.csv',row.names=F)
