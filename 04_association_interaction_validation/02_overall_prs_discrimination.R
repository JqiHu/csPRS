# Analyze prediction accuracy of PRSs
setwd('/path/to/csPRS/score_split_analysis/analysis')
library(data.table)
library(survival)
library(stringr)
library(survcomp)

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

# Cindex
getCindex <- function(p){
  prs.name.tmp <- paste0(p,c('_csbeta','_beta'))
  # C-index
  res <- c()
  for(exp in prs.name.tmp){
    cindex <- concordance.index(x=dat[,exp],surv.time=dat$CAD.inci.fu,
				surv.event=dat$CAD.inci)
    add <- paste0(round(cindex$c.index,2),
                 ' (',round(cindex$lower,2),
                 '-',round(cindex$upper,2),')')

    add <- c(exp,add)
    res <- rbind(res,add)
  }

  # Clean output
  res <- cbind(p,res)

  return(res)
  print(paste0('Finished for ',p))
}

#### Run across PRSs
# clusters
cluster <- str_remove_all(prs_name[grep('_csbeta',prs_name)],'_csbeta')
res.list <- lapply(cluster,getCindex)
res <- do.call(rbind,res.list)
colnames(res)[1:2] <- c('cluster','type')

write.csv(res,'cindex_inciCAD_compare.csv',row.names=F)
