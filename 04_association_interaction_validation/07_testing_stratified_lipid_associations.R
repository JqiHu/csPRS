# Stratified association between lipids and CAD by subgroup membership
setwd('/path/to/csPRS/score_split_analysis/analysis_test/subgroup_tertile')
library(data.table)
library(survival)
library(stringr)

# Step 1. Prepare data
## phenotypes
dat.phe <- fread('../../data/ukb_inciCAD.txt.gz',data.table=F)
## subgroup information
subgroup <- fread('subgroup.tsv',data.table=F)

dat <- merge(dat.phe,subgroup,by='eid',all.x=T)
###### !!!!!!!!!!!!!!!!!!!!!
###### Analyze among testing samples!!!!!!!!!
dat <- dat[dat$batch=='testing',]
dim(dat)

# Add lipids information
lipids <- fread('/path/to/resources/phenotypes/data/ukb_lipids_related.csv',data.table=F)
##### organize data
colnames(lipids) <- c('eid','med1','med2','ApolipoproteinA','ApolipoproteinB',
                      'cholesterol','HDL','LDL','LipoproteinA','Triglycerides')
#### subjects without med to lower cholesterol
row.keep <- union(which(lipids$med1 != 1),which(lipids$med2 != 1))
#lipids.add <- lipids[row.keep,]
lipids.add <- lipids
lipids.add[-row.keep,c('ApolipoproteinA','ApolipoproteinB','cholesterol','HDL','LDL','LipoproteinA','Triglycerides')] <- NA
lipids.add$med <- 0
lipids.add$med[which(lipids$med1==1|lipids$med2==1)] <- 1
(table(lipids.add$med))

dat <- merge(dat,lipids.add,by='eid')
print(dim(dat))
print(table(dat$CAD.inci))

##### Define binary lipids based on normal ranges
## Convention: 0 = within normal range, 1 = abnormal
## Apolipoprotein A (ApoA) – normal: ≥ 1 g/L
dat$apoa <- ifelse(dat$ApolipoproteinA >= 1.00, 0, 1)
print(table(dat$apoa, useNA = "ifany"))
## Apolipoprotein B (ApoB) – normal:
## Males: ≤ 1.33 g/L; Females: ≤ 1.17 g/L
dat$apob <- with(dat,
                 ifelse((sex == 1 & ApolipoproteinB <= 1.33) |
                        (sex == 0 & ApolipoproteinB <= 1.17),
                        0, 1))
print(table(dat$apob, useNA = "ifany"))
## Total cholesterol – normal: < 5.17 mmol/L
dat$chol <- ifelse(dat$cholesterol < 5.17, 0, 1)
print(table(dat$chol, useNA = "ifany"))
## HDL – low HDL is abnormal
## Males: normal if ≥ 1.03 mmol/L
## Females: normal if ≥ 1.29 mmol/L
dat$hdl_bin <- with(dat,
                    ifelse((sex == 1 & HDL < 1.03) |
                           (sex == 0 & HDL < 1.29),
                           1, 0))
print(table(dat$hdl_bin, useNA = "ifany"))
## LDL-C – normal: < 2.59 mmol/L
dat$ldl_bin <- ifelse(dat$LDL < 2.59, 0, 1)
print(table(dat$ldl_bin, useNA = "ifany"))
## Lipoprotein(a) – normal: < 75 nmol/L
dat$lipa <- ifelse(dat$LipoproteinA < 75, 0, 1)
print(table(dat$lipa, useNA = "ifany"))
## Triglycerides – normal: < 1.69 mmol/L
dat$tg <- ifelse(dat$Triglycerides < 1.69, 0, 1)
print(table(dat$tg, useNA = "ifany"))
lipids.names.ref <- data.frame(quant=colnames(lipids.add)[4:11],
                           bi=c('apoa','apob','chol','hdl_bin',
                                'ldl_bin','lipa','tg','med'))

# Step 2. Stratified associations between lipids and CAD
#### Function of associations
getHR <- function(g,l){
  # extract variables
  covar <- c('sex','age_recruit',paste0('PC',1:10))
  tmp <- subset(dat,select=c('CAD.inci',g,l,covar,'CAD.inci.fu'))
  colnames(tmp)[2:3] <- c('subgroup','lipids')

  ### Create three outcome
  ## 1. among High csPRS
  tmp$out1[tmp$subgroup=='High' & tmp$CAD.inci==0] <- 0
  tmp$out1[tmp$subgroup=='High' & tmp$CAD.inci==1] <- 1 
  ## 2. among intermediate csPRS
  tmp$out2[tmp$subgroup=='Intermediate' & tmp$CAD.inci==0] <- 0
  tmp$out2[tmp$subgroup=='Intermediate' & tmp$CAD.inci==1] <- 1
  ## 3. among low csPRS
  tmp$out3[tmp$subgroup=='Low' & tmp$CAD.inci==0] <- 0
  tmp$out3[tmp$subgroup=='Low' & tmp$CAD.inci==1] <- 1
  print(table(tmp$out2,tmp$out1,useNA='ifany'))
  print(table(tmp$out3,tmp$out1,useNA='ifany'))

  ### associations stratified 
  formu <- formula(paste0('Surv(CAD.inci.fu,CAD.inci) ~ lipids+',paste(covar,collapse='+')))
  ##### outcome1
  fit1 <- coxph(formu,data=tmp[is.na(tmp$out1)==F,])
  coef1 <- summary(fit1)$coefficients
  add1 <- coef1['lipids',]

  ##### outcome2
  fit2 <- coxph(formu,data=tmp[is.na(tmp$out2)==F,])
  coef2 <- summary(fit2)$coefficients
  add2 <- coef2['lipids',]

  ##### outcome3
  fit3 <- coxph(formu,data=tmp[is.na(tmp$out3)==F,])
  coef3 <- summary(fit3)$coefficients
  add3 <- coef3['lipids',]

  res <- rbind(add1,add2,add3)
  res <- cbind(g,l,c('High','Intermediate','Low'),res)

  return(res)
}

# Run for all PRSs and lipids
lipids.names <- lipids.names.ref$bi
group <- colnames(subgroup)[-c(1,2)]
pairs <- expand.grid(group, lipids.names, stringsAsFactors = FALSE)
res <- mapply(getHR, pairs$Var1, pairs$Var2, SIMPLIFY = FALSE)
res <- do.call(rbind, res)

res <- as.data.frame(res)
colnames(res)[1:3] <- c('subgroup','lipids','type')
res$lipids_ori <- lipids.names.ref$quant[match(res$lipids,lipids.names.ref$bi)]
res$subgroup <- str_remove_all(res$subgroup,'_csbeta')
res$coef <- as.numeric(res$coef)
res$`se(coef)` <- as.numeric(res$`se(coef)`)

res$HR <- exp(res$coef)
res$HR_low <- exp(res$coef - 1.96*res$`se(coef)`)
res$HR_up <- exp(res$coef + 1.96*res$`se(coef)`)

write.csv(res,'subgroup_lipids_stratify.csv',row.names=F)


