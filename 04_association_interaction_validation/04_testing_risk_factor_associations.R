# Associations with CAD risk factors
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

# Add lipids information
lipids <- fread('/path/to/resources/phenotypes/data/ukb_lipids_related.csv',data.table=F)
##### organize data
colnames(lipids) <- c('eid','med1','med2','ApolipoproteinA','ApolipoproteinB',
                      'cholesterol','HDL','LDL','LipoproteinA','Triglycerides')
#### subjects without med to lower cholesterol
row.keep <- union(which(lipids$med1 != 1),which(lipids$med2 != 1)) # no overlap between med1 and med2
#lipids.add <- lipids[row.keep,]
lipids.add <- lipids
#### For subjects who are taking medicine
##### set the lipids level to missing!!!!!!!!!
lipids.add[-row.keep,c('ApolipoproteinA','ApolipoproteinB','cholesterol','HDL','LDL','LipoproteinA','Triglycerides')] <- NA
## add use of lipids-lowering medicine
lipids.add$med <- 0
lipids.add$med[which(lipids$med1==1|lipids$med2==1)] <- 1
(table(lipids.add$med))

dat <- merge(dat,lipids.add,by='eid')

# Risk factors for CAD
dat.exp <- fread('/path/to/resources/phenotypes/data_clean/ukb_CAD_risk_factor.tsv.gz',data.table=F)
#### extract non-lipids
lipid <- c('HDL','LDL','TC','Apolipoprotein_A',
           'Apolipoprotein_B','Lipoprotein_A','logTG')
phe.extract <- setdiff(colnames(dat.exp),lipid)
dat.exp2 <- subset(dat.exp,select=phe.extract)

# Merge data
covar0 <- c('age_recruit','sex',paste0('PC',1:10))
dat <- merge(dat,dat.exp2,by='eid')
phe.all <- colnames(dat)[45:115]

# Function to get OR for specific PD-PRSs
getOR <- function(p){
  ## cell type
  cell <- str_remove_all(p,'_csbeta')
  covar0 <- c('age_recruit','sex',paste0('PC',1:10))

  res <- c()
  for(o in phe.all){
    tmp <- subset(dat,select=c(o,p,covar0))
    colnames(tmp)[1:2] <- c('outcome','prs')
    tmp <- tmp[complete.cases(tmp),]
    if(o %in% c('BrC','OvC')){
      tmp <- tmp[tmp$sex==0,]
      covar <- setdiff(covar0,'sex')
    }else if(o == 'PrC'){
      tmp <- tmp[tmp$sex==1,]
      covar <- setdiff(covar0,'sex')
    }else{
      covar <- covar0
    }

    formu <- formula(paste0('outcome ~ prs+',paste(covar,collapse='+')))
    ### check variable type
    if(length(unique(tmp$outcome))==2){
      fit <- glm(formu,data=tmp,family='binomial')
    }else{
      fit <- lm(formu,data=tmp)
    }
    coef <- summary(fit)$coefficients
    add <- coef['prs',]
    add <- c(o,add)

    res <- rbind(res,add)
  }
  res <- cbind.data.frame(cell,res)
  colnames(res)[3] <- 'outcome'
  print(paste0('Finished for ',cell))
  return(res)
}

#### Run across PRSs
# clusters
res.list <- lapply(prs_name,getOR)
res <- do.call(rbind,res.list)
colnames(res)[1] <- 'PRS'

write.csv(res,'test_or_risk_factor.csv',row.names=F)
