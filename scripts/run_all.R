# ============================================================
# ARDS cell-state signature — reproducibility script
# Run from repository root: Rscript scripts/run_all.R
# ============================================================
setwd("/data/wenjiying/ARDS/ARDS_repo")
set.seed(123)

suppressPackageStartupMessages({
  library(pROC)
  library(GSVA)
  library(dplyr)
})

# ============================================================
# Step 1: Build sample ledger
# ============================================================
cat("=== Step 1: Build sample ledger ===\n")
ledger <- read.csv("/data/wenjiying/ARDS/analysis/GSE32707_ARDS_vs_Untreated_ledger.csv")
cat("  Total:", nrow(ledger),
    " ARDS:", sum(ledger$disease == "ARDS"),
    " Untreated:", sum(ledger$disease == "Untreated"), "\n")
write.csv(ledger, "data_processed/ledger.csv", row.names = FALSE)

# ============================================================
# Step 2: Load ssGSEA scores and define labels
# ============================================================
cat("\n=== Step 2: Load ssGSEA scores ===\n")
scores <- as.matrix(read.csv("/data/wenjiying/ARDS/analysis/ssGSEA_32707.csv",
                             row.names = 1, check.names = FALSE))
common <- intersect(colnames(scores), ledger$geo_accession)
scores <- scores[, common, drop = FALSE]
labels <- ifelse(ledger$disease[match(common, ledger$geo_accession)] == "ARDS", 1, 0)
key <- c("Immature_neutrophils", "Non-classical_monocytes")
cat("  Samples:", length(labels), " ARDS:", sum(labels), "\n")

# ============================================================
# Step 3: Nested cross-validation
# ============================================================
cat("\n=== Step 3: Nested cross-validation ===\n")
df <- data.frame(t(scores[key, , drop = FALSE]), label = labels)
colnames(df) <- c(key, "label")
cat("  df dim:", dim(df), " colnames:", colnames(df), "\n")

n_rep <- 10
n_fold <- 5
all_preds <- c()
all_labels <- c()
for (r in 1:n_rep) {
  idx <- sample(1:nrow(df))
  fid <- integer(nrow(df))
  fid[idx] <- cut(seq_along(idx), n_fold, labels = FALSE)
  for (f in 1:n_fold) {
    te <- which(fid == f)
    tr <- which(fid != f)
    fit <- glm(label ~ ., data = df[tr, , drop = FALSE], family = binomial)
    all_preds <- c(all_preds, predict(fit, newdata = df[te, , drop = FALSE], type = "response"))
    all_labels <- c(all_labels, df$label[te])
  }
}
roc_cv <- roc(all_labels, all_preds, quiet = TRUE)
cv_auc <- as.numeric(auc(roc_cv))
ci <- ci.auc(roc_cv)
cat("  CV AUC:", round(cv_auc, 3),
    " (95% CI", round(ci[1], 3), "-", round(ci[3], 3), ")\n")
write.csv(data.frame(AUC = cv_auc, CI_low = ci[1], CI_high = ci[3]),
          "results/Fig2_CV_AUC.csv", row.names = FALSE)

# ============================================================
# Step 4: Train fixed model
# ============================================================
cat("\n=== Step 4: Train fixed model ===\n")
fit_final <- glm(label ~ ., data = df, family = binomial)
coef_df <- data.frame(feature = names(coef(fit_final)),
                      coefficient = as.numeric(coef(fit_final)))
write.csv(coef_df, "model/model_coefficients.csv", row.names = FALSE)
saveRDS(fit_final, "model/model_fixed.rds")
cat("  Coefficients:\n")
print(coef_df)

pred_all <- predict(fit_final, type = "response")
write.csv(data.frame(sample = common,
                     group = ifelse(labels == 1, "ARDS", "Untreated"),
                     predicted_probability = pred_all),
          "results/per_sample_predictions.csv", row.names = FALSE)

# ============================================================
# Step 5: External validation
# ============================================================
cat("\n=== Step 5: External validation ===\n")
validate_one <- function(scores_mat, labels_vec, name) {
  common_key <- intersect(key, rownames(scores_mat))
  df_test <- data.frame(t(scores_mat[common_key, , drop = FALSE]))
  colnames(df_test) <- common_key
  pred <- predict(fit_final, newdata = df_test, type = "response")
  roc_obj <- roc(labels_vec, pred, quiet = TRUE)
  ci <- ci.auc(roc_obj)
  cat("  ", name, "AUC:", round(auc(roc_obj), 3),
      " (", round(ci[1], 3), "-", round(ci[3], 3), ")\n")
  data.frame(cohort = name, AUC = as.numeric(auc(roc_obj)),
             CI_low = ci[1], CI_high = ci[3])
}

res_list <- list()

s243 <- as.matrix(read.csv("/data/wenjiying/ARDS/analysis/ssGSEA_243066.csv",
                           row.names = 1, check.names = FALSE))
l243 <- ifelse(grepl("^A", colnames(s243)), 1, 0)
res_list[[1]] <- validate_one(s243, l243, "GSE243066")

s200 <- as.matrix(read.csv("/data/wenjiying/ARDS/analysis/ssGSEA_GSE200847.csv",
                           row.names = 1, check.names = FALSE))
gm <- c(TA203="Hypo",TA213="Hypo",TA225="Hypo",TA227="Hyper",TA229="Hyper",
  TA234="Hypo",TA235="Hypo",TA251="Hypo",TA257="Hypo",TA288="Hypo",
  TA311="Hyper",TA314="Hyper",TA315="Hypo",TA319="Control",TA320="Hyper",
  TA334="Hypo",TA335="Hyper",TA337="Hyper",TA341="Hypo",TA353="Hyper",
  TA358="Hypo",TA392="Hypo",TA400="Hypo",TA402="Hypo",TA409="Hypo",
  TA411="Hypo",TA413="Hypo",TA414="Hypo",TA458="Hypo",TA462="Control",
  TA463="Hypo",TA477="Control",TA479="Hypo",TA488="Hypo",TA503="Hyper",
  TA524="Hypo",TA530="Hypo",TA539="Hypo",TA540="Hypo",TA541="Hypo",
  TA543="Hypo",TA551="Hypo",TA555="Control",TA563="Control",TA565="Hyper",
  TA575="Hypo")
l200 <- ifelse(gm[colnames(s200)] %in% c("Hyper", "Hypo"), 1, 0)
res_list[[2]] <- validate_one(s200, l200, "GSE200847")

s171 <- as.matrix(read.csv("/data/wenjiying/ARDS/analysis/pb171524_scores.csv",
                           row.names = 1, check.names = FALSE))
l171 <- ifelse(grepl("^L", colnames(s171)), 1, 0)
res_list[[3]] <- validate_one(s171, l171, "GSE171524")

res_ext <- do.call(rbind, res_list)
write.csv(res_ext, "results/Fig3_external_AUC.csv", row.names = FALSE)

# ============================================================
# Step 6: Threshold optimization
# ============================================================
cat("\n=== Step 6: Threshold optimization ===\n")
roc_obj <- roc(labels, pred_all, quiet = TRUE)
best <- coords(roc_obj, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
ac <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
ac <- ac[ac$threshold > 0 & ac$threshold < 1, ]
sens95 <- ac[ac$sensitivity >= 0.95, ]
if (nrow(sens95) > 0) sens95 <- sens95[which.max(sens95$specificity), ]
spec90 <- ac[ac$specificity >= 0.90, ]
if (nrow(spec90) > 0) spec90 <- spec90[which.max(spec90$sensitivity), ]
out <- rbind(
  data.frame(type = "Youden", best[, c("threshold", "sensitivity", "specificity")]),
  if (exists("sens95") && nrow(sens95) > 0)
    data.frame(type = "HighSens", sens95[, c("threshold", "sensitivity", "specificity")]),
  if (exists("spec90") && nrow(spec90) > 0)
    data.frame(type = "HighSpec", spec90[, c("threshold", "sensitivity", "specificity")])
)
write.csv(out, "results/Fig6_clinical_utility.csv", row.names = FALSE)
cat("  Thresholds:\n")
print(out)

# ============================================================
# Step 7: Spatial statistics
# ============================================================
cat("\n=== Step 7: Spatial statistics ===\n")
sp <- read.csv("/data/wenjiying/ARDS/analysis/spatial_summary.csv")
write.csv(sp, "results/Fig4_spatial_stats.csv", row.names = FALSE)
cat("  Spatial summary:\n")
print(sp)

cat("\n=== ALL DONE ===\n")
