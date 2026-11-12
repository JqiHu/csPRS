# Associations with incident CAD
# Among testing samples
setwd('/path/to/csPRS/score_split_analysis/analysis_test')
library(data.table)
library(survival)
library(stringr)

# Step 1. Prepare data
## phenotypes
dat.phe <- fread('../data/ukb_inciCAD.txt.gz',data.table=F)
## PRSs
dat.prs <- fread('../data/ukb_clustercsprs.txt.gz',data.table=F)
### focuse on csPRSs
dat.prs <- subset(dat.prs,select=c('eid',colnames(dat.prs)[grep('csbeta',colnames(dat.prs))]))
prs_name <- colnames(dat.prs)[-1]
### genome-wide PRS
dat.prs1 <- fread('/path/to/resources/prs/score/SDPR_cad_cardiogramc4a.profile',
                  data.table=F,select=c('IID','SCORE'))
colnames(dat.prs1) <- c('eid','PRS_CAD')

dat <- merge(dat.phe,dat.prs,by='eid')
dat <- merge(dat,dat.prs1,by='eid')
prs_name <- c(prs_name,'PRS_CAD')
###### !!!!!!!!!!!!!!!!!!!!!
###### Analyze among testing samples!!!!!!!!!
dat <- dat[dat$batch=='testing',]
dim(dat)
### standardize PRSs
dat[,prs_name] <- apply(dat[,prs_name],2,function(x){as.numeric(scale(x))})

### Hazards ratios of each PRS of the same cell type on CAD
getHR <- function(p){
  covar <- c('age_recruit','sex',paste0('PC',1:10))

  # fit cox
  formu <- formula(paste0('Surv(CAD.inci.fu,CAD.inci) ~ ',p,'+',paste(covar,collapse='+')))
  fit <- coxph(formu,data=dat)
  coef <- summary(fit)$coefficients
  add <- coef[p,]

  res <- c(p,add)
  print(paste0('Finished for ',p))
  return(res)
}

#### Run across PRSs
# clusters
res.list <- lapply(prs_name,getHR)
res <- do.call(rbind,res.list)
colnames(res)[1] <- 'PRS'

res <- as.data.frame(res)
res$PRS <- str_remove_all(res$PRS,'_csbeta')
res$coef <- as.numeric(res$coef)
res$`se(coef)` <- as.numeric(res$`se(coef)`)
res$HR_low <- exp(res$coef-1.96*res$`se(coef)`)
res$HR_up <- exp(res$coef+1.96*res$`se(coef)`)

write.csv(res,'test_hr_inciCAD.csv',row.names=F)
