# Traing XGBoost model
setwd('/path/to/csPRS/score_split_analysis/combine/xgboost')
library(data.table)
library(survival)
library(xgboost)

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

# Step 3. Prepare survival labels ----------
# One common encoding for 'survival:cox' in xgboost:
#   label = time, with negative sign for censored observations
label_train <- ifelse(train$CAD.inci == 1, train$CAD.inci.fu, -train$CAD.inci.fu)

dtrain <- xgb.DMatrix(data = x_train_scaled, label = label_train)
dtest  <- xgb.DMatrix(data = x_test_scaled)

# Step 4. Specify parameters & cross-validated training ----------
params <- list(
  objective        = "survival:cox",
  eval_metric      = "cox-nloglik",

  # --- Tree complexity ---
  max_depth        = 2,        # shallow trees → avoid overfitting
  min_child_weight = 1,        # small because data likely not huge
  gamma            = 0,        # regularization for splits

  # --- Learning behavior ---
  eta              = 0.05,     # low learning rate
  subsample        = 0.8,      # row subsampling
  colsample_bytree = 0.8,      # feature subsampling
  colsample_bylevel= 0.8,

  # --- Regularization ---
  lambda           = 1,
  alpha            = 0
)

set.seed(116)
cv_xgb <- xgb.cv(
  params                = params,
  data                  = dtrain,
  nrounds               = 500,
  nfold                 = 5,
  early_stopping_rounds = 20,
  verbose               = 0
)

best_nrounds <- cv_xgb$best_iteration

xgb_fit <- xgb.train(
  params  = params,
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

# Step 5. Predictions & C-index on TEST ----------
xgb_risk_test <- as.numeric(predict(xgb_fit, newdata = dtest))

c_xgb <- concordance(y_test ~ xgb_risk_test)$concordance

# Step 6. Feature importance ----------
importance_mat <- xgb.importance(
  model = xgb_fit,
  feature_names = colnames(x_train_scaled)
)
# importance_mat has Gain, Cover, Frequency, etc.

# Step 7. Save everything ----------
# 7.1 Model object
saveRDS(xgb_fit, file = "xgb_cox_model.rds")

# 7.2 Test predictions with outcomes & IDs
xgb_pred_test <- data.frame(
  id         = test[['eid']],
  time       = test$CAD.inci.fu,
  status     = test$CAD.inci,
  xgb_risk   = xgb_risk_test
)
write.csv(xgb_pred_test, "xgb_cox_test_predictions.csv", row.names = FALSE)

# 7.3 Feature importance
write.csv(importance_mat, "xgb_cox_importance.csv", row.names = FALSE)

# 7.4 Metrics summary
metrics_xgb <- data.frame(
  model      = "xgb_cox",
  c_index    = c_xgb,
  nrounds    = best_nrounds,
  n_train    = nrow(train),
  n_test     = nrow(test)
)
write.csv(metrics_xgb, "xgb_cox_metrics.csv", row.names = FALSE)

