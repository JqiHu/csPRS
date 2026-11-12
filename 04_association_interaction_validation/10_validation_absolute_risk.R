# Absolute interactions
setwd('/path/to/csPRS/score_split_analysis/analysis_vali/subgroup_tertile')
library(data.table)
library(stringr)
library(dplyr)
library(survival)

# Step 1. Prepare data
## phenotypes
dat.phe <- fread('../../data/ukb_inciCAD.txt.gz',data.table=F)
## subgroup information
subgroup <- fread('subgroup.tsv',data.table=F)

dat <- merge(dat.phe,subgroup[,-2],by='eid')
###### !!!!!!!!!!!!!!!!!!!!!
dat <- dat[dat$batch=='vali',]
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

lipids.names <- data.frame(quant=colnames(lipids.add)[4:11],
			   bi=c('apoa','apob','chol','hdl_bin',
                  	        'ldl_bin','lipa','tg','med'))

#------------------------------------------
# Step 2. Absolute interactions
#------------------------------------------
# Function to get AR, ARR, and test for ARR trend
getAR <- function(l,p){
  tmp <- subset(dat,select=c(l,p,'CAD.inci','CAD.inci.fu'))
  colnames(tmp)[1:2] <- c('l','p')
  tmp <- tmp[is.na(tmp$l)==F,]

  ## Define six groups
  tmp$group <- paste0(tmp$p,'_',tmp$l)
  print(table(tmp$group))

  ## AR in each group
  ## AR=proportion of CAD
  tab <- tmp %>%
    group_by(group) %>%
    summarise(
      n  = n(),
      events = sum(CAD.inci),
      IR = events / n,
      .groups = "drop"
    )

  ## ---- add IRR by strata ----
  tab$IRR <- NA
  tab$IRR[tab$group=='High_1'] <- tab$IR[tab$group=='High_1']-tab$IR[tab$group=='High_0']
  tab$IRR[tab$group=='Intermediate_1'] <- tab$IR[tab$group=='Intermediate_1']-tab$IR[tab$group=='Intermediate_0']
  tab$IRR[tab$group=='Low_1'] <- tab$IR[tab$group=='Low_1']-tab$IR[tab$group=='Low_0']
   
  ## ---- final table ----
  results <- cbind.data.frame(l,p,tab)

  ## Test for trend
  tmp$p <- factor(tmp$p,levels=c('Low','Intermediate','High'))
  tmp$p_num <- as.numeric(tmp$p)
  fit <- coxph(Surv(CAD.inci.fu, CAD.inci) ~ l * p_num,
               data = tmp)
  p.trend <- summary(fit)$coefficients['l:p_num',5]

  results$p.trend <- p.trend

  return(results)
}

# Run for PRS and lipids of interest
VAR1 <- colnames(subgroup)[-c(1:2)]
VAR2 <- lipids.names$bi
pairs <- expand.grid(VAR1 = VAR1,
                     VAR2 = VAR2,
                     stringsAsFactors = FALSE)
res <- mapply(getAR, pairs$VAR2, pairs$VAR1, SIMPLIFY = FALSE)
res <- do.call(rbind, res)

res <- as.data.frame(res)
res$subgroup <- str_remove_all(res$p,'_csbeta')
res$lipi <- lipids.names$quant[match(res$l,lipids.names$bi)]

write.csv(res,'subgroup_lipids_inter_absolute.csv',row.names=F)


