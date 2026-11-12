# Associations with incident CAD
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

### Hazards ratios of each PRS of the same cell type on CAD
getHR <- function(p){
  covar <- c('age_recruit','sex',paste0('PC',1:10))
  prs.names.tmp <- paste0(p,c('_csbeta','_beta'))

  # fit cox
  res <- c()
  for(exp in prs.names.tmp){
    formu <- formula(paste0('Surv(CAD.inci.fu,CAD.inci) ~ ',exp,'+',paste(covar,collapse='+')))
    fit <- coxph(formu,data=dat)
    coef <- summary(fit)$coefficients
    add <- coef[exp,]

    res <- rbind(res,add)
  }

  # Clean output
  res <- cbind(p,prs.names.tmp,res)

  return(res)
  print(paste0('Finished for ',p))
}

#### Run across PRSs
# clusters
cluster <- str_remove_all(prs_name[grep('_csbeta',prs_name)],'_csbeta')
res.list <- lapply(cluster,getHR)
res <- do.call(rbind,res.list)
colnames(res)[1:2] <- c('cluster','type')

res <- as.data.frame(res)
res$coef <- as.numeric(res$coef)
res$`se(coef)` <- as.numeric(res$`se(coef)`)
res$HR_low <- exp(res$coef-1.96*res$`se(coef)`)
res$HR_up <- exp(res$coef+1.96*res$`se(coef)`)

write.csv(res,'hr_inciCAD_compare.csv',row.names=F)
