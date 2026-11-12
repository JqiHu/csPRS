# Interactions with lipids for the incidence of CAD
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

# Function of interactions
getInter <- function(prs,lipid){
  covar <- c('age_recruit','sex',paste0('PC',1:10))
  # extract data
  tmp <- subset(dat,select=c('eid','CAD.inci','CAD.inci.fu',prs,lipid,covar))
  colnames(tmp)[4:5] <- c('prs','lipid')
  tmp <- tmp[complete.cases(tmp),]
  ### create categorical PRS
  # Get tertile cut points
  cuts <- quantile(tmp$prs, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  # Create tertile variable
  tmp$prs_cate <- cut(
    tmp$prs,
    breaks = cuts,
    include.lowest = TRUE,
    labels = c("T1", "T2", "T3")
  )
  tmp$prs_cate <- factor(tmp$prs_cate,levels=c('T1','T2','T3'),
			 label=c('lowPRS','midPRS','highPRS'))

  # fit cox
  # Full model (with interaction)
  formu_full <- formula(
    paste0(
      'Surv(CAD.inci.fu, CAD.inci) ~ prs_cate*lipid + ',
      paste(covar, collapse = '+')
    )
  )

  fit_full <- coxph(formu_full, data = tmp)

  # Reduced model (no interaction)
  formu_red <- formula(
    paste0(
      'Surv(CAD.inci.fu, CAD.inci) ~ prs_cate + lipid + ',
      paste(covar, collapse = '+')
    )
  )

  fit_red <- coxph(formu_red, data = tmp)

  # Likelihood ratio test
  lrt <- anova(fit_red, fit_full, test = "LRT")

  # Extract p-value
  p_lrt <- lrt[2, "Pr(>|Chi|)"]

  ### Extract interaction terms
  coef_table <- summary(fit_full)$coefficients
  conf_table <- summary(fit_full)$conf.int

  # Keep only interaction rows
  interaction_rows <- grep(":", rownames(coef_table))

  out <- data.frame(
    prs = prs,
    lipids = lipid,
    Term = rownames(coef_table)[interaction_rows],
    HR   = conf_table[interaction_rows, "exp(coef)"],
    Lower95 = conf_table[interaction_rows, "lower .95"],
    Upper95 = conf_table[interaction_rows, "upper .95"],
    Wald_p  = coef_table[interaction_rows, "Pr(>|z|)"],
    row.names = NULL
  )
  out$lrt_p <- p_lrt
  return(out)

}

# Run for all PRSs and lipids
lipids.names <- lipids.names.ref$bi
pairs <- expand.grid(prs_name, lipids.names, stringsAsFactors = FALSE)
res <- mapply(getInter, pairs$Var1, pairs$Var2, SIMPLIFY = FALSE)
res <- do.call(rbind, res)

res <- as.data.frame(res)
res$lipids.ori <- lipids.names.ref$quant[match(res$lipids,lipids.names.ref$bi)]

write.csv(res,'test_inter_lipids_prs_cate_lipi_cate.csv',row.names=T)
