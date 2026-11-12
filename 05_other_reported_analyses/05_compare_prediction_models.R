# Compare C-index of CAD PRS and the integrated PRSs
setwd('/path/to/csPRS/score_split_analysis/combine/res')
library(data.table)
library(survival)
library(dplyr)
library(survcomp)

# Step 1. Prepare data
## phenotypes
dat <- fread('../../data/ukb_inciCAD.txt.gz',data.table=F)
test  <- dat[dat$batch=='testing', ]

### genome-wide PRS
dat.prs1 <- fread('/path/to/resources/prs/score/SDPR_cad_cardiogramc4a.profile',
                  data.table=F,select=c('IID','SCORE'))
colnames(dat.prs1) <- c('eid','PRS_CAD')

test <- merge(test,dat.prs1,by='eid')

### Read predictions from ridge and xgboost
ridge <- fread('../ridge/ridge_cox_test_predictions.csv',data.table=F,
	       select=c('id','ridge_lp'))
enet <- fread('../elastic_net/enet_cox_test_predictions.csv',data.table=F,
	      select=c('id','ridge_lp'))
colnames(enet)[2] <- 'enet_lp'
xgb <- fread('../xgboost/xgb_cox_test_predictions.csv',data.table=F,
	     select=c('id','xgb_risk'))

test <- merge(test,ridge,by.x='eid',by.y='id')
test <- merge(test,enet,by.x='eid',by.y='id')
test <- merge(test,xgb,by.x='eid',by.y='id')

# Cindex
out <- c()
for(p in c('PRS_CAD','ridge_lp','enet_lp','xgb_risk')){
  cindex <- concordance.index(x=test[,p],surv.time=test$CAD.inci.fu,
			      surv.event=test$CAD.inci)
  add <- c(p,cindex$c.index,cindex$lower,cindex$upper)
  
  out <- rbind(out,add)
}

out <- as.data.frame(out)
colnames(out) <- c('model','cindex','lower','upper')
out[,-1] <- apply(out[,-1],2,as.numeric)
out$value <- paste0(round(out$cindex,2),' (',
		    round(out$lower,2),'-',
		    round(out$upper,2),')')

write.csv(out,'cindex_compare.csv')
