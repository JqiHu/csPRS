# Combine csPRSs using Ridge regression 
setwd('/path/to/csPRS/score_split_analysis/combine/ridge')
library(data.table)
library(survival)
library(glmnet)

# Step 1. Prepare data
## phenotypes
dat.phe <- fread('../../data/ukb_inciCAD.txt.gz',data.table=F)
## PRSs
### csPRSs
dat.prs <- fread('../../data/ukb_clustercsprs.txt.gz',data.table=F)
### remaining csPRSs
dat.prs2 <- fread('../../data/ukb_clustercsprs_remain.txt.gz',data.table=F)

###### Focus only on csBETA
dat.prs <- subset(dat.prs,select=c('eid',colnames(dat.prs)[grep('csbeta',colnames(dat.prs))]))
dat.prs2 <- subset(dat.prs2,select=c('eid',colnames(dat.prs2)[grep('csbeta',colnames(dat.prs2))]))
colnames(dat.prs2)[-1] <- paste0('remain_',colnames(dat.prs2)[-1])

dat.prs <- merge(dat.prs,dat.prs2,by='eid')

prs_name <- colnames(dat.prs)[-1]

dat <- merge(dat.phe,dat.prs,by='eid')

train <- dat[dat$batch=='training', ]
test  <- dat[dat$batch=='testing', ]

y_train <- Surv(train$CAD.inci.fu, train$CAD.inci)
y_test  <- Surv(test$CAD.inci.fu,  test$CAD.inci)

# Step 2. Build design matrices (scale using TRAIN only) 
x_train_raw <- as.matrix(train[, prs_name])
x_test_raw  <- as.matrix(test[, prs_name])

### Training standardization
x_train_scaled <- scale(x_train_raw)
train_mean  <- attr(x_train_scaled, "scaled:center")
train_sd    <- attr(x_train_scaled, "scaled:scale")

### Testing standardization — use its own mean/sd
x_test_scaled <- scale(x_test_raw)
test_mean  <- attr(x_test_scaled, "scaled:center")
test_sd    <- attr(x_test_scaled, "scaled:scale")

# Step 3. Fit ridge Cox with cross-validation ----------
set.seed(116)
cv_ridge <- cv.glmnet(
  x      = x_train_scaled,
  y      = y_train,
  family = "cox",
  alpha  = 0       # ridge
)

best_lambda <- cv_ridge$lambda.min  # or lambda.1se

# Step 4. Predictions & performance on TEST ----------
ridge_lp_test <- as.numeric(predict(
  cv_ridge,
  newx = x_test_scaled,
  s    = "lambda.min",
  type = "link"       # linear predictor (risk score)
))

c_ridge <- concordance(y_test ~ ridge_lp_test)$concordance

# Step 5. Extract feature importance (coefficients) ----------
coef_mat <- as.matrix(coef(cv_ridge, s = "lambda.min"))
coef_df  <- data.frame(
  feature = rownames(coef_mat),
  coef    = as.numeric(coef_mat)
)
rownames(coef_df) <- NULL

# Step 6. Save everything for later ----------
# 6.1 Model object
saveRDS(cv_ridge, file = "ridge_cox_model.rds")

# 6.2 Test predictions with outcomes and IDs
ridge_pred_test <- data.frame(
  id        = test[['eid']],
  time      = test$CAD.inci.fu,
  status    = test$CAD.inci,
  ridge_lp  = ridge_lp_test
)
write.csv(ridge_pred_test, "ridge_cox_test_predictions.csv", row.names = FALSE)

# 6.3 Coefficients
write.csv(coef_df, "ridge_cox_coefficients.csv", row.names = FALSE)

# 6.4 Metrics summary
metrics_ridge <- data.frame(
  model    = "ridge_cox",
  c_index  = c_ridge,
  lambda   = best_lambda,
  n_train  = nrow(train),
  n_test   = nrow(test)
)
write.csv(metrics_ridge, "ridge_cox_metrics.csv", row.names = FALSE)
