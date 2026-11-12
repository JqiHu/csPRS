# Correlations across csPRSs
setwd('/path/to/csPRS/score_split_analysis/analysis')
library(data.table)
library(survival)
library(stringr)

# Step 1. Prepare data
## phenotypes
dat.phe <- fread('../data/ukb_inciCAD.txt.gz',data.table=F)
## PRSs
dat.prs <- fread('../data/ukb_clustercsprs.txt.gz',data.table=F)
prs_name <- colnames(dat.prs)[-1]

dat <- merge(dat.phe,dat.prs,by='eid')
###### !!!!!!!!!!!!!!!!!!!!!
###### Analyze among training samples!!!!!!!!!
dat <- dat[dat$batch=='training',]
dim(dat)
### standardize PRSs
dat[,prs_name] <- apply(dat[,prs_name],2,function(x){as.numeric(scale(x))})

library(Hmisc)
library(dplyr)
prs <- dat[,prs_name]
### pair-wise correlations
res <- rcorr(as.matrix(prs), type = "pearson")
# Correlation matrix
cor_matrix <- res$r
# P-value matrix
pval_matrix <- res$P

# Convert to long format
cor_long <- melt(res$r, varnames = c("Var1", "Var2"), value.name = "correlation")
pval_long <- melt(res$P, varnames = c("Var1", "Var2"), value.name = "p_value")

# Merge correlation and p-value tables
long_df <- left_join(cor_long, pval_long, by = c("Var1", "Var2"))

write.csv(long_df,'prs_corr.csv',row.names=F)



