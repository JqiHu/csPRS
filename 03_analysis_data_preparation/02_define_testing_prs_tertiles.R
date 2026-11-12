# Define subgroups among patients 
# while compare cases with subjects with similar genetic risk
setwd('/path/to/csPRS/score_split_analysis/analysis_test/subgroup_tertile')
library(data.table)

# Step 1. Prepare data
## phenotypes
dat.phe <- fread('../../data/ukb_inciCAD.txt.gz',data.table=F)
## PRSs
dat.prs <- fread('../../data/ukb_clustercsprs.txt.gz',data.table=F)
### focuse on csPRSs
dat.prs <- subset(dat.prs,select=c('eid',colnames(dat.prs)[grep('csbeta',colnames(dat.prs))]))
prs_name <- colnames(dat.prs)[-1]
### genome-wide PRS
dat.prs1 <- fread('/path/to/resources/prs/score/SDPR_cad_cardiogramc4a.profile',
                  data.table=F,select=c('IID','SCORE'))
colnames(dat.prs1) <- c('eid','PRS_CAD')

dat <- merge(dat.phe,dat.prs,by='eid')
dat <- merge(dat,dat.prs1,by='eid')
#prs_name <- c(prs_name,'PRS_CAD')
###### !!!!!!!!!!!!!!!!!!!!!
###### Analyze among testing samples!!!!!!!!!
dat2 <- dat[dat$batch=='testing',]
dim(dat2)

# Function to get tertiles
getTer <- function(score){
  # Get tertile cut points
  cuts <- quantile(dat2[[score]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  # Create tertile variable
  group <- cut(dat2[[score]],breaks=cuts,include.lowest=T,
	       	labels=c('Low','Intermediate','High'))
  out <- data.frame(eid=dat2$eid,group=group)
  colnames(out)[2] <- score
  return(out)
}

#### run for each csPRS
out <- data.frame(eid=dat2$eid,CAD=dat2$CAD.inci)
prs_name <- c('PRS_CAD',prs_name)
for(p in prs_name){
  # get subgroup subjects ID
  subgroup <- getTer(p)

  out <- merge(out,subgroup,by='eid',all.x=T)
  colnames(out)[ncol(out)] <- p
}

summary(out)
write.table(out,'subgroup.tsv',
            row.names=F,sep='\t',quote=F)

