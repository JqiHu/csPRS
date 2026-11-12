# Quantify disease trajectory using embeddings
setwd('/path/to/csPRS/score_split_analysis/analysis_test/subgroup_tertile')
library(reticulate)
library(data.table)

# read subgroup definitions
subgroup <- fread('subgroup.tsv',data.table=F)
### keep CAES only
subgroup <- subgroup[subgroup$CAD==1,]

# covariates
covar0 <- c('age_recruit','sex',paste0('PC',1:10))
dat.cov <- fread('/path/to/resources/phenotypes/data/ukbb_kinship_eur_PCs.csv',
                 select=c('eid',covar0),data.table=F)
subgroup <- merge(dat.cov,subgroup,by='eid')

# Embeddings supplied by collaborators
pickle   <- import("pickle")
builtins <- import("builtins")

path <- "/path/to/resources/embeddings/disease_embeddings.pickle"

f <- builtins$open(path, "rb")
obj <- pickle$load(f)
f$close()
obj.icd10 <- names(obj)

# Risk factors for CAD with TIME
dat.exp <- fread('/path/to/resources/phenotypes/data_clean/ukb_CAD_risk_factor_wTime.tsv.gz',data.table=F)
## start with diseases
ref <- list('AF'=c(paste0('I48',c(0:4,9))),
	    'Angina'=paste0('I20',c(0,1,8,9)),
	    'HF'=c('I110','I130','I132','I255','I420','I425',
                   'I428','I429','I500','I1501','I1509','I150'),
	    'IS'=c(paste0('I60',0:9),
                   paste0('I61',c(0:6,8,9)),
                   paste0('I62',c(0,1,9)),
                   paste0('I63',c(0:6,8,9)),
                   'I64'),
	    'PAD'=c('I7000','I7001','I7020','I7021','I7080','I7081','I7090','I7091',
                'I738','I739'),
	    'PAH'=c('I270'),
	    'T2D'=c(paste0('E11',0:9)),
	    'CKD'=c('D593','D638','E102','E112','E132','E142',
               paste0('E85',c(0:4,8,9)),
               paste0('I12',c(8,9)),
               paste0('I13',c(1,2)),
               'I151','M103','N001','N011','N020','N021',
               'N030','N031',paste0('N04',c(0:9)),
               'N050','N051','N060','N061','N071',
               paste0('N08',c(0:5)),
               'N10',paste0('N11',c(0,1,8,9)),'N12',
               paste0('N14',c(1:4)),'N165',
               paste0('N17',c(0:2,8,9)),
               paste0('N18',c(0:5,9)),'N19',
               'N250','N990',paste0('Q61',c(0:5,8,9)),
               'T824','T861','Y602','Y841',paste0('Z49',c(0:2)),
               'Z940','Z992')
	    )
ref.icd10 <- unlist(ref)
diseases <- names(ref)
setdiff(ref.icd10,obj.icd10) # "I1501" "I1509" "I128"  "N046" 

dat.exp2 <- subset(dat.exp,select=c('eid','CAD','CAD_days',diseases,
                                    paste0(diseases,'_days')))
# Merge
dat <- merge(subgroup,dat.exp2,by='eid')
dat <- as.data.table(dat)

#############
## 1. Keep only diseases diagnosed before CAD
#############
cad_col <- "CAD_days"
dat <- as.data.table(dat)

# CAD time
stopifnot(cad_col %in% names(dat))
dat[, (cad_col) := as.numeric(get(cad_col))]
dat[get(cad_col) <= 0, (cad_col) := NA_real_]
dat <- dat[!is.na(get(cad_col))]   # should keep nearly all

# disease columns
d_cols <- paste0(diseases, "_days")
missing_cols <- setdiff(d_cols, names(dat))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse=", "))

for (cc in d_cols) {
  dat[, (cc) := as.numeric(get(cc))]
  dat[get(cc) < 0, (cc) := NA_real_]
}
dat[, `:=`(T = get(cad_col), event = 1L)]  # time-to-CAD, event indicator

## Clean each disease relative to CAD (pre-CAD window)
for (d in diseases) {
  td <- paste0(d, "_days")

  # event before CAD
  dat[, (paste0("event_", d, "_preCAD")) :=
        fifelse(!is.na(get(td)) & get(td) < T, 1L, 0L)]

  # time within pre-CAD window:
  # if diagnosed before CAD => actual time
  # else => set to CAD time (meaning "not yet by CAD")
  dat[, (paste0("time_", d, "_preCAD")) :=
        fifelse(!is.na(get(td)) & get(td) < T, get(td), T)]
}

#############
## 2. Individual-embedding
#############
# helper: mean-embedding across a set of ICD10 codes
mean_embed_codes <- function(codes, obj) {
  codes <- unique(codes)
  codes <- codes[codes %in% names(obj)]   # keep only codes that exist in obj
  if (length(codes) == 0) return(NULL)

  mats <- lapply(codes, function(cd) obj[[cd]])
  # convert to matrix: rows=codes, cols=embedding dims
  M <- do.call(rbind, mats)
  colMeans(M)
}

# compute disease embeddings (named list)
disease_embed_list <- lapply(names(ref), function(d) mean_embed_codes(ref[[d]], obj))
names(disease_embed_list) <- names(ref)

# drop diseases with no available ICD10 embeddings
keep <- !vapply(disease_embed_list, is.null, logical(1))
disease_embed_list <- disease_embed_list[keep]

# stack into a matrix: diseases x dims
disease_mat <- do.call(rbind, disease_embed_list)
# rownames(disease_mat) are diseases

# --------------------------------------
#### Individual embedding as average
# find the event columns and map them to disease names
event_cols <- grep("^event_.*_preCAD$", names(dat), value = TRUE)
event_diseases <- sub("^event_(.*)_preCAD$", "\\1", event_cols)

# keep only diseases that have embeddings
common <- intersect(event_diseases, rownames(disease_mat))
stopifnot(length(common) > 0)

# subset dat + disease_mat to the same disease order
dat <- as.data.frame(dat)
event_cols2 <- paste0("event_", common, "_preCAD")
eid <- dat$eid
E <- as.matrix(dat[, c(event_cols2)])

# coerce logical -> numeric if needed
if (is.logical(E)) E <- 1L * E

# weighted sum of disease embeddings for each subject: (n x diseases) %*% (diseases x dims) = n x dims
S <- E %*% disease_mat[common, , drop = FALSE]

# divide by number of pre-CAD diseases for each subject to get the mean
k <- rowSums(E)
indiv_mat <- S / pmax(k, 1)  # avoid divide-by-0

# if no pre-CAD diseases, set embedding to NA
indiv_mat[k == 0, ] <- NA_real_

# output with eid
out <- data.frame(
  eid = eid,
  as.data.frame(indiv_mat, check.names = FALSE),
  check.names = FALSE
)
colnames(out)[-1] <- paste0("emb_", seq_len(ncol(indiv_mat)))
sum(is.na(out$emb_1)) # 3,284 missing...

# indicator: has any pre-CAD disease
#has_predisease <- as.integer(k > 0)
# fill embeddings for k==0 with 0 instead of NA
#indiv_mat_filled <- indiv_mat
#indiv_mat_filled[k == 0, ] <- NA

# output with eid + embeddings + binary dimension
#out <- data.frame(
#  eid = eid,
#  as.data.frame(indiv_mat_filled, check.names = FALSE),
#  check.names = FALSE
#)
#colnames(out)[2:(ncol(indiv_mat_filled) + 1)] <- paste0("emb_", seq_len(ncol(indiv_mat_filled)))


#############
## 3. Compare embedding across subgroups
#############
subgroup_names <- colnames(dat)[grep('csbeta',colnames(dat))]
dat2 <- merge(out,subgroup,by='eid')

# keep only subjects with embeddings
emb_cols <- grep("^emb_", names(dat2), value = TRUE)
dat2 <- dat2[complete.cases(dat2[, emb_cols]), ]
dim(dat2) # 1,495 subjects
apply(dat2[,subgroup_names],2,table)

# PERMANOVA
library(vegan)

X <- as.matrix(dat2[, c(emb_cols)])
D <- dist(X)
##### decrease dimensions by PCA
#pc <- prcomp(X, scale. = TRUE)

#var <- cumsum(pc$sdev^2 / sum(pc$sdev^2))
#k <- which(var >= 0.8)[1]   # e.g., 80%
#k

#X_red <- pc$x[, 1:k]
#D_red <- dist(X_red)

covar <- c('age_recruit','sex',paste0('PC',1:10))
res_list <- lapply(subgroup_names, function(s) {

  # build formula
  formu <- as.formula(
    paste("D ~", paste(c(covar,s), collapse = " + "))
  )

  fit <- adonis2(
    formu,
    data = dat2,
    permutations = 1000
  )

  out <- as.data.frame(fit)
  out$term <- rownames(out)

  # extract subgroup row only
  sg <- out[out$term == s, , drop = FALSE]

  data.frame(
    subgroup = s,
    R2 = sg$R2,
    p_value = sg$`Pr(>F)`
  )
})

# combine results
perm_res <- do.call(rbind, res_list)
write.csv(perm_res,'embedding_caseOnly.csv',row.names=F)
