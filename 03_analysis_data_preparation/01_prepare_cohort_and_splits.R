# Prepare incident CAD
# and splitting
library(data.table)
library(stringr)

# Phenotype
dat <- fread('/path/to/resources/phenotypes/data/ukbb_CAD.tsv',
             data.table=F)
# relatedness and PCs from UKB
cov <- fread('/path/to/resources/phenotypes/data/ukbb_kinship_eur_PCs.csv',
             data.table=F,select=c('eid','kin','eur',paste0('PC',1:10)))
# predicted ancestry supplied by a collaborator
ancestry <- fread('/path/to/resources/ancestry/ukb_ancestry.tsv',
                  data.table=F,select=c('IID','predicted_pop'))

dat <- merge(dat,cov,by='eid')
dat <- merge(dat,ancestry,by.x='eid',by.y='IID')
print(nrow(dat)) # 487164

# --------------------------------------------------
### First, keep only unrelated subjecst
dat.tmp1 <- dat[dat$kin==0,]
print(nrow(dat.tmp1)) # 339490
table(dat.tmp1$predicted_pop,dat.tmp1$eur,useNA='always')

# --------------------------------------------------
### Second, remove all non-EUR
dat.tmp2 <- dat.tmp1[dat.tmp1$predicted_pop=='EUR',]
print(nrow(dat.tmp2)) # 318217
table(dat.tmp2$predicted_pop,dat.tmp2$eur,useNA='always')

# --------------------------------------------------
### Third, remove prevalent CAD cases
dat.w <- dat.tmp2
###### Create incident cases!!!!!!!
dat.w$CAD.inci[dat.w$CAD==0] <- 0
dat.w$CAD.inci[dat.w$CAD==1] <- 1
table(dat.w$CAD,dat.w$CAD.inci,useNA='always')
###### follow-up time from recruit to end
dat.w$CAD.inci.fu <- as.Date(dat.w$CAD_date_end)-as.Date(dat.w$date_recruit)
dat.w$CAD.inci.fu <- as.numeric(dat.w$CAD.inci.fu)
range(dat.w$CAD.inci.fu[is.na(dat.w$CAD.inci)==F]) # missing!
### check missing
sum(is.na(dat.w$CAD.inci.fu[is.na(dat.w$CAD.inci)==F]))
missing <- dat.w[is.na(dat.w$CAD.inci.fu)&is.na(dat.w$CAD.inci)==F,]
######## impute end date
# Split into year and month
parts <- strsplit(missing$CAD_date_end, "\\.")
# Convert to data.frame for easier handling
df <- do.call(rbind, lapply(parts, function(p) {
  year <- as.integer(p[1])
  month <- ifelse(length(p) > 1, as.integer(p[2]), 1L)  # default to Jan if missing
  data.frame(year, month)
}))
# Build Date with day = 1
missing$CAD_date_end2 <- as.Date(sprintf("%04d-%02d-01", df$year, df$month))
missing$CAD.inci.fu <- as.numeric(missing$CAD_date_end2 - as.Date(missing$date_recruit))

dat.w$CAD.inci.fu[match(missing$eid,dat.w$eid)] <- missing$CAD.inci.fu

#### Some negative values
sum(dat.w$CAD.inci.fu[is.na(dat.w$CAD.inci)==F]<0)
dat.w$CAD.inci[dat.w$CAD.inci.fu<0] <- NA
table(dat.w$CAD,dat.w$CAD.inci,useNA='always')
range(dat.w$CAD.inci.fu[is.na(dat.w$CAD.inci)==F])
###### Keep only incident CAD
dat.w <- dat.w[is.na(dat.w$CAD.inci)==F,]
print(nrow(dat.w)) # 308072
print(table(dat.w$CAD.inci,dat.w$eur,useNA='always'))

# --------------------------------------------------
#### Last, divide dataset into training (60% Caucasian), testing (40% Caucasian), and validation (non-Cauc EUR)
dat <- dat.w
dat$batch[is.na(dat$eur)] <- 'vali'

set.seed(116)
n <- sum(is.na(dat$eur)==F)

idx <- sample(n, size = 0.6 * n)  # 60% indices
dat$batch[is.na(dat$eur)==F][idx] <- 'training'
dat$batch[is.na(dat$batch)] <- 'testing'

print(table(dat$batch,dat$eur,useNA='ifany'))

### Check number of CAD 
print(table(dat$CAD.inci,dat$batch))

setwd('/path/to/csPRS/score_split_analysis/data/')
write.table(dat,gzfile('ukb_inciCAD.txt.gz'),
	    row.names=F,sep='\t',quote=F)


