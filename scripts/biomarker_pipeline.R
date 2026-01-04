# ─────────────────────────────────────────────────────────────────────────────
# Biomarker Discovery Pipeline 
# Train: GSE9006 (GPL96, whole blood)  →  Validate: GSE15932 (GPL570, blood)
# Outputs:
#   - discovery_consensus_ranking.csv
#   - discovery_panel_genes.txt
#   - panel_coefficients.csv
#   - table_external_performance.csv
#   - table_confusion_matrix_counts.csv
#   - table_confusion_matrix_rates.csv
#   - table_panel_genes.csv
#   - calibration_validation.png
#   - figure1_external_roc_ggplot.png/pdf
#   - pca_pre_vs_post_combined.png
#   - biomarker_discovery_outputs.rds
# ─────────────────────────────────────────────────────────────────────────────

## 0) Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("tidyverse","caret","pROC","randomForest","glmnet","Boruta","ggplot2","patchwork")
for (p in cran) if (!require(p, character.only = TRUE)) install.packages(p)
bio <- c("GEOquery","limma","hgu133a.db","hgu133plus2.db","AnnotationDbi")
for (p in bio) if (!require(p, character.only = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)

set.seed(123)

## 1) Helper functions
ensure_log2 <- function(mat) {
  qs <- quantile(mat, probs = c(.25,.5,.75), na.rm = TRUE)
  if (isTRUE(qs[3] > 100)) mat <- log2(mat + 1)
  mat
}
youden_threshold <- function(roc_obj) {
  coords(roc_obj, "best", best.method = "youden",
         ret = c("threshold","sensitivity","specificity"))
}
ppv_npv <- function(sens, spec, prev) {
  ppv <- (sens*prev) / (sens*prev + (1-spec)*(1-prev))
  npv <- (spec*(1-prev)) / ((1-sens)*prev + spec*(1-prev))
  c(PPV = ppv, NPV = npv)
}

# ─────────────────────────────────────────────────────────────────────────────
# 2) TRAIN: GSE9006 (GPL96) → gene-symbol matrix
# ─────────────────────────────────────────────────────────────────────────────
message("Downloading GSE9006 (training)...")
gse  <- GEOquery::getGEO("GSE9006", GSEMatrix = TRUE)[[1]]
expr_probe <- exprs(gse)
expr_probe <- ensure_log2(expr_probe)

meta <- pData(gse)
meta$Status <- ifelse(grepl("type 2", meta$characteristics_ch1.4, ignore.case = TRUE),
                      "T2D","Control")
metadata_train <- meta |>
  dplyr::transmute(SampleID = geo_accession,
                   Status   = factor(Status, levels = c("Control","T2D"))) |>
  dplyr::filter(!is.na(Status))

# Map probes→SYMBOL on GPL96 (hgu133a)
probe2sym <- AnnotationDbi::mapIds(hgu133a.db,
                                   keys = rownames(expr_probe),
                                   column = "SYMBOL",
                                   keytype = "PROBEID",
                                   multiVals = "first")

expr_sym <- as.data.frame(expr_probe) |>
  tibble::rownames_to_column("ProbeID") |>
  dplyr::mutate(Gene = unname(probe2sym[ProbeID])) |>
  dplyr::filter(!is.na(Gene)) |>
  dplyr::select(-ProbeID) |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(dplyr::across(everything(), mean), .groups = "drop") |>
  as.data.frame()
rownames(expr_sym) <- expr_sym$Gene; expr_sym$Gene <- NULL

# Align columns to metadata
common_train <- intersect(colnames(expr_sym), metadata_train$SampleID)
expr_sym <- expr_sym[, common_train, drop = FALSE]
metadata_train <- metadata_train[match(common_train, metadata_train$SampleID), , drop = FALSE]

# ─────────────────────────────────────────────────────────────────────────────
# 3) Nested CV for DEGs and fold-stability (limma within folds)
# ─────────────────────────────────────────────────────────────────────────────
message("Running fold-stability DEG selection (nested CV)...")
outer_folds <- caret::createFolds(metadata_train$Status, k = 5, returnTrain = TRUE)
deg_bag   <- list()
auc_nested <- numeric(length(outer_folds))

for (i in seq_along(outer_folds)) {
  tr_ids   <- outer_folds[[i]]
  te_ids   <- setdiff(seq_len(nrow(metadata_train)), tr_ids)
  tr_samps <- metadata_train$SampleID[tr_ids]
  te_samps <- metadata_train$SampleID[te_ids]
  
  expr_tr <- expr_sym[, tr_samps, drop = FALSE]
  grp     <- factor(metadata_train$Status[tr_ids], levels = c("Control","T2D"))
  design  <- model.matrix(~ grp)
  
  fit  <- limma::eBayes(limma::lmFit(expr_tr, design))
  degs <- limma::topTable(fit, coef = 2, number = Inf, adjust = "BH")
  keep <- subset(degs, adj.P.Val < 0.05 & abs(logFC) > 1)
  topg <- head(rownames(keep), 50)  # up to 50 per fold
  deg_bag[[i]] <- keep[topg, c("logFC","adj.P.Val")]
  
  # Quick nested ROC using RF on fold-top genes (if enough)
  if (length(topg) >= 5) {
    tr_df <- data.frame(t(expr_tr[topg,,drop=FALSE]), Status = grp)
    te_df <- data.frame(t(expr_sym[topg, te_samps, drop=FALSE]),
                        Status = metadata_train$Status[te_ids])
    ctrl  <- caret::trainControl(method = "cv", number = 3, classProbs = TRUE,
                                 summaryFunction = caret::twoClassSummary)
    rf_m  <- caret::train(Status ~ ., data = tr_df,
                          method = "rf", trControl = ctrl, metric = "ROC")
    pr    <- predict(rf_m, newdata = te_df, type = "prob")[,"T2D"]
    auc_nested[i] <- pROC::auc(pROC::roc(te_df$Status, pr))
  } else {
    auc_nested[i] <- NA_real_
  }
}
cat("Nested CV mean AUC (RF on fold-top genes):",
    round(mean(auc_nested, na.rm = TRUE), 3), "\n")

# Stability & effect size summary
stab_tbl <- purrr::imap_dfr(deg_bag, ~{
  tibble::tibble(Gene = rownames(.x), logFC = .x$logFC)
})
stab_summary <- stab_tbl |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(sel_freq = n() / length(outer_folds),
                   mean_abs_logFC = mean(abs(logFC)), .groups = "drop")

# ─────────────────────────────────────────────────────────────────────────────
# 4) Training frame on union of stable genes
# ─────────────────────────────────────────────────────────────────────────────
stable_pool <- stab_summary |>
  dplyr::arrange(dplyr::desc(sel_freq), dplyr::desc(mean_abs_logFC)) |>
  dplyr::pull(Gene) |>
  unique()

if (length(stable_pool) < 20) {
  message("Few stable DEGs; falling back to top variable genes.")
  vr <- apply(expr_sym, 1, sd)
  stable_pool <- names(sort(vr, decreasing = TRUE))[1:500]
}

X_tr <- t(expr_sym[stable_pool, , drop = FALSE])
y_tr <- metadata_train$Status
df_tr <- data.frame(X_tr, Status = y_tr)
y_bin <- ifelse(df_tr$Status == "T2D", 1, 0)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Model-agnostic importance — Boruta + bootstrapped elastic-net
# ─────────────────────────────────────────────────────────────────────────────
message("Running Boruta...")
bor <- Boruta::Boruta(Status ~ ., data = df_tr, maxRuns = 200, doTrace = 0)
bor_imp <- Boruta::TentativeRoughFix(bor)
bor_stats <- Boruta::attStats(bor_imp)
bor_calls <- data.frame(Gene = rownames(bor_stats),
                        boruta_call = bor_stats$decision, row.names = NULL)

message("Bootstrapped elastic-net selection (B = 500)...")
B <- 500
enet_sel <- vector("list", B)
x_mat <- as.matrix(df_tr[, colnames(df_tr) != "Status"])

for (b in 1:B) {
  idx <- sample(seq_len(nrow(x_mat)), replace = TRUE)
  x_b <- x_mat[idx, , drop = FALSE]
  y_b <- y_bin[idx]
  cv  <- glmnet::cv.glmnet(x_b, y_b, family = "binomial",
                           alpha = 0.5, nfolds = 5, type.measure = "auc")
  fit <- glmnet::glmnet(x_b, y_b, family = "binomial",
                        alpha = 0.5, lambda = cv$lambda.1se)
  genes_b <- rownames(as.matrix(fit$beta))[as.numeric(fit$beta) != 0]
  enet_sel[[b]] <- genes_b
}
enet_prob <- enet_sel |>
  unlist() |>
  table() |>
  as.data.frame() |>
  setNames(c("Gene","count")) |>
  dplyr::mutate(prob = count / B) |>
  dplyr::select(Gene, prob)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Consensus ranking & small interpretable panel
# ─────────────────────────────────────────────────────────────────────────────
consensus <- stab_summary |>
  dplyr::left_join(bor_calls, by = "Gene") |>
  dplyr::left_join(enet_prob, by = "Gene") |>
  dplyr::mutate(
    prob = dplyr::coalesce(prob, 0),
    boruta_flag = ifelse(boruta_call == "Confirmed", 1, 0),
    score = 0.5*sel_freq + 0.4*prob + 0.1*as.numeric(scale(mean_abs_logFC))
  ) |>
  dplyr::arrange(dplyr::desc(score))

panel <- consensus |>
  dplyr::filter(sel_freq >= 0.6, (prob >= 0.6 | boruta_flag == 1)) |>
  dplyr::arrange(dplyr::desc(score)) |>
  dplyr::slice_head(n = 20) |>
  dplyr::pull(Gene)

if (length(panel) < 5) {
  panel <- consensus$Gene[1:min(15, nrow(consensus))]
}
cat("Chosen panel size:", length(panel), "\n")
getwd()
write.csv(consensus, "discovery_consensus_ranking1.csv", row.names = FALSE)
writeLines(panel, "discovery_panel_genes1.txt")

# ─────────────────────────────────────────────────────────────────────────────
# 7) Fit a shrunken, interpretable model (elastic-net/logistic)
# ─────────────────────────────────────────────────────────────────────────────
x_panel <- as.matrix(df_tr[, panel, drop = FALSE])
cv_pan  <- glmnet::cv.glmnet(x_panel, y_bin, family = "binomial",
                             alpha = 0.5, nfolds = 5, type.measure = "auc")
fit_pan <- glmnet::glmnet(x_panel, y_bin, family = "binomial",
                          alpha = 0.5, lambda = cv_pan$lambda.1se)

cf <- as.matrix(coef(fit_pan))
coefs <- data.frame(Term = rownames(cf), Coef = as.numeric(cf[,1]), row.names = NULL)
write.csv(coefs, "panel_coefficients.csv", row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 8) EXTERNAL VALIDATION: GSE15932 (GPL570)
# ─────────────────────────────────────────────────────────────────────────────
message("Downloading GSE15932 (validation)...")
eset_val <- GEOquery::getGEO("GSE15932", GSEMatrix = TRUE)[[1]]
expr_val_probe <- ensure_log2(exprs(eset_val))
pd_val <- pData(eset_val)
tt <- tolower(pd_val$title)

# keep only T2D and healthy, remove pancreatic cancer
is_cancer <- grepl("pancreatic cancer", tt)
is_t2d    <- grepl("patient\\s*b\\d+", tt) | grepl("diabetes mellitus", tt)
is_ctrl   <- grepl("healthy control", tt)
keep_idx  <- (is_t2d | is_ctrl) & !is_cancer

pd_val        <- pd_val[keep_idx, , drop = FALSE]
expr_val_probe<- expr_val_probe[, rownames(pd_val), drop = FALSE]
Status_val <- ifelse(grepl("healthy control", tolower(pd_val$title)), "Control", "T2D")
y_va <- factor(Status_val, levels = c("Control","T2D"))

# Map probes→SYMBOL on GPL570 (hgu133plus2)
probe2sym570 <- AnnotationDbi::mapIds(hgu133plus2.db,
                                      keys = rownames(expr_val_probe),
                                      column = "SYMBOL", keytype = "PROBEID",
                                      multiVals = "first")
expr_sym_val <- as.data.frame(expr_val_probe) |>
  tibble::rownames_to_column("ProbeID") |>
  dplyr::mutate(Gene = unname(probe2sym570[ProbeID])) |>
  dplyr::filter(!is.na(Gene)) |>
  dplyr::select(-ProbeID) |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(dplyr::across(everything(), mean), .groups = "drop") |>
  as.data.frame()
rownames(expr_sym_val) <- expr_sym_val$Gene; expr_sym_val$Gene <- NULL
expr_sym_val <- as.matrix(expr_sym_val)

# Intersection with panel
common_panel <- intersect(panel, rownames(expr_sym_val))
if (length(common_panel) < 5) stop("Too few panel genes present in GSE15932.")
X_tr_pan <- as.matrix(df_tr[, common_panel, drop = FALSE])
X_va_pan <- t(expr_sym_val[common_panel, , drop = FALSE])

# Batch-correct (train vs val) while preserving class
X_merge <- rbind(X_tr_pan, X_va_pan)
Status_all <- factor(c(as.character(y_tr), as.character(y_va)), levels = c("Control","T2D"))
batch <- factor(c(rep("train", nrow(X_tr_pan)), rep("val", nrow(X_va_pan))))
X_adj <- t(limma::removeBatchEffect(t(X_merge), batch = batch,
                                    design = model.matrix(~ Status_all)))

X_tr_adj <- X_adj[seq_len(nrow(X_tr_pan)), , drop = FALSE]
X_va_adj <- X_adj[-seq_len(nrow(X_tr_pan)), , drop = FALSE]

# Refit elastic-net on adjusted train; predict validation
cv_pan2 <- glmnet::cv.glmnet(X_tr_adj, y_bin, family = "binomial",
                             alpha = 0.5, nfolds = 5, type.measure = "auc")
fit_pan2 <- glmnet::glmnet(X_tr_adj, y_bin, family = "binomial",
                           alpha = 0.5, lambda = cv_pan2$lambda.1se)
prob_va  <- as.numeric(predict(fit_pan2, newx = X_va_adj, type = "response"))

roc_va   <- pROC::roc(y_va, prob_va)
auc_va   <- pROC::auc(roc_va)
ci_va    <- pROC::ci.auc(roc_va, boot.n = 2000)
cat(sprintf("External AUC (panel) on GSE15932: %.3f | 95%% CI: %.3f–%.3f\n",
            auc_va, ci_va[1], ci_va[3]))

# ─────────────────────────────────────────────────────────────────────────────
# 9) Thresholded clinical metrics (Youden), PPV/NPV at assumed prevalence
# ─────────────────────────────────────────────────────────────────────────────
thr <- youden_threshold(roc_va)
thr$youden <- with(thr, sensitivity + specificity - 1)
best_idx <- order(thr$youden, thr$sensitivity, thr$specificity, decreasing = TRUE)[1]
best <- thr[best_idx, , drop = FALSE]

prev_assumed <- 0.11
p  <- as.numeric(prev_assumed)
Se <- as.numeric(best$sensitivity)
Sp <- as.numeric(best$specificity)

PPV <- (Se * p) / (Se * p + (1 - Sp) * (1 - p))
NPV <- (Sp * (1 - p)) / (Sp * (1 - p) + (1 - Se) * p)

cat(sprintf(
  "Chosen (Youden-max, sensitivity-priority): threshold=%.6f | Sens=%.3f | Spec=%.3f | Youden=%.3f\n",
  best$threshold, Se, Sp, Se + Sp - 1
))
cat(sprintf("At prevalence %.0f%% → PPV=%.3f | NPV=%.3f\n", p * 100, PPV, NPV))

# ─────────────────────────────────────────────────────────────────────────────
# 10) Calibration & Brier + plot
# ─────────────────────────────────────────────────────────────────────────────
brier <- mean((as.numeric(y_va == "T2D") - prob_va)^2)
cat(sprintf("Brier score (validation) = %.3f\n", brier))

cal_df <- data.frame(y = as.numeric(y_va == "T2D"), p = prob_va) |>
  dplyr::mutate(bin = cut(p, breaks = quantile(p, probs = seq(0,1,0.1)),
                          include.lowest = TRUE)) |>
  dplyr::group_by(bin) |>
  dplyr::summarise(pred = mean(p), obs = mean(y), n = dplyr::n(), .groups = "drop")

ggplot(cal_df, aes(x = pred, y = obs, size = n)) +
  geom_point() + geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(title = "Calibration (validation)", x = "Predicted risk", y = "Observed T2D rate") +
  theme_minimal()
ggsave("calibration_validation.png", width = 5, height = 4, dpi = 150)

# ─────────────────────────────────────────────────────────────────────────────
# 11) Permutation test (label shuffling on validation)
# ─────────────────────────────────────────────────────────────────────────────
message("Permutation test (label shuffling on validation)...")
K <- 500
perm_auc <- numeric(K)
for (k in 1:K) {
  y_perm <- sample(y_va)
  perm_auc[k] <- as.numeric(pROC::auc(pROC::roc(y_perm, prob_va)))
}
p_perm <- (1 + sum(perm_auc >= as.numeric(auc_va))) / (1 + K)
cat(sprintf("Permutation p (AUC vs chance) = %.4f\n", p_perm))

# ─────────────────────────────────────────────────────────────────────────────
# 12) Save ROC figure (ggplot) with stats annotations
# ─────────────────────────────────────────────────────────────────────────────
auc_val <- as.numeric(auc_va)
ci_auc  <- ci_va
p_text  <- sprintf("AUC = %.3f (95%% CI %.3f–%.3f)\nPermutation p = %.4f (K = 500)",
                   auc_val, ci_auc[1], ci_auc[3], p_perm)

g <- pROC::ggroc(roc_va, legacy.axes = TRUE) +
  labs(title = "External ROC — GSE15932 (Diabetes vs Control)",
       x = "False Positive Rate", y = "True Positive Rate") +
  annotate("text", x = 0.65, y = 0.15, label = p_text, size = 3.5, hjust = 0) +
  theme_minimal()

ggsave("figure1_external_roc_ggplot.pdf", g, width = 6.5, height = 5)
ggsave("figure1_external_roc_ggplot.png", g, width = 6.5, height = 5, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# 13) Performance tables & confusion matrix at chosen threshold
# ─────────────────────────────────────────────────────────────────────────────
pred_class <- factor(ifelse(prob_va >= as.numeric(best$threshold), "T2D", "Control"),
                     levels = c("Control","T2D"))
tab <- table(Observed = y_va, Predicted = pred_class)
TN <- tab["Control","Control"]; FP <- tab["Control","T2D"]
FN <- tab["T2D","Control"];     TP <- tab["T2D","T2D"]

sens   <- ifelse((TP+FN)>0, TP/(TP+FN), NA_real_)
spec   <- ifelse((TN+FP)>0, TN/(TN+FP), NA_real_)
youden <- sens + spec - 1
acc    <- (TP+TN)/sum(tab)
balacc <- mean(c(sens, spec), na.rm = TRUE)

# CIs for Sens/Spec (Wilson via prop.test)
sens_ci <- suppressWarnings(prop.test(TP, TP+FN)$conf.int)
spec_ci <- suppressWarnings(prop.test(TN, TN+FP)$conf.int)

perf_tbl <- data.frame(
  Metric = c("AUC","AUC_CI_low","AUC_CI_high","Brier",
             "Youden_threshold","Sensitivity","Specificity",
             "Sensitivity_CI_low","Sensitivity_CI_high",
             "Specificity_CI_low","Specificity_CI_high",
             "Youden_index","Accuracy","Balanced_accuracy",
             "PPV_at_prev","NPV_at_prev","Prevalence_assumed",
             "Permutation_p"),
  Value  = c(as.numeric(auc_va), as.numeric(ci_va[1]), as.numeric(ci_va[3]), brier,
             as.numeric(best$threshold), sens, spec,
             sens_ci[1], sens_ci[2], spec_ci[1], spec_ci[2],
             youden, acc, balacc, PPV, NPV, p, p_perm)
)
write.csv(perf_tbl, "table_external_performance.csv", row.names = FALSE)

cm_counts <- as.data.frame.matrix(tab)
cm_rates  <- data.frame(
  Metric = c("Sensitivity (TPR)","Specificity (TNR)","PPV","NPV","Accuracy","Balanced accuracy"),
  Value  = c(sens, spec, PPV, NPV, acc, balacc)
)
write.csv(cm_counts, "table_confusion_matrix_counts.csv")
write.csv(cm_rates,  "table_confusion_matrix_rates.csv", row.names = FALSE)

panel_tbl <- consensus |>
  dplyr::filter(Gene %in% panel) |>
  dplyr::left_join(
    coefs |>
      dplyr::rename(Gene = Term, Coefficient = Coef),
    by = "Gene"
  ) |>
  dplyr::filter(Gene != "(Intercept)") |>
  dplyr::transmute(
    Gene,
    Selection_Frequency = round(sel_freq, 2),
    ENet_Selection_Prob = round(dplyr::coalesce(prob, 0), 2),
    Boruta = dplyr::if_else(boruta_flag == 1, "Confirmed", "Rejected/Tentative"),
    Mean_abs_logFC = round(mean_abs_logFC, 2),
    Consensus_Score = round(score, 3),
    Coefficient = round(Coefficient, 3),
    Direction_vs_T2D = ifelse(Coefficient >= 0, "Higher in T2D", "Lower in T2D")
  ) |>
  dplyr::arrange(dplyr::desc(Consensus_Score))
write.csv(panel_tbl, "table_panel_genes.csv", row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 14) PCA before vs after batch correction (combined figure)
# ─────────────────────────────────────────────────────────────────────────────
pca_plot <- function(mat, batch, status, tag){
  pc  <- prcomp(mat, center = TRUE, scale. = TRUE)
  var <- (pc$sdev^2) / sum(pc$sdev^2) * 100
  df  <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2],
                    Batch = factor(batch, levels = c("train","val")),
                    Status = factor(status, levels = c("Control","T2D")))
  ggplot(df, aes(PC1, PC2, colour = Batch, shape = Status)) +
    geom_point(size = 2.6, alpha = 0.9) +
    labs(title = paste("PCA", tag, "batch correction"),
         x = sprintf("PC1 (%.1f%%)", var[1]),
         y = sprintf("PC2 (%.1f%%)", var[2])) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}
p_pre  <- pca_plot(X_merge, batch, Status_all, "before")
p_post <- pca_plot(X_adj,   batch, Status_all, "after")

combined <- p_pre | p_post
combined <- combined + plot_annotation(
  title = "PCA before vs after batch correction",
  subtitle = "Color = Batch (train/val), Shape = Status (Control/T2D)"
)
ggsave("pca_pre_vs_post_combined.png", combined, width = 11, height = 5.5, dpi = 150)

# ─────────────────────────────────────────────────────────────────────────────
# 15) Save all key R objects to RDS
# ─────────────────────────────────────────────────────────────────────────────
saveRDS(list(
  panel = panel, consensus = consensus, coefs = coefs,
  roc = roc_va, thr = thr, brier = brier, perm_p = p_perm,
  auc_nested = auc_nested, X_tr_adj = X_tr_adj, X_va_adj = X_va_adj,
  Status_all = Status_all, batch = batch
), file = "biomarker_discovery_outputs.rds")

message("Done. Key files:\n",
        "  • discovery_consensus_ranking.csv\n",
        "  • discovery_panel_genes.txt\n",
        "  • panel_coefficients.csv\n",
        "  • table_external_performance.csv\n",
        "  • table_confusion_matrix_counts.csv\n",
        "  • table_confusion_matrix_rates.csv\n",
        "  • table_panel_genes.csv\n",
        "  • calibration_validation.png\n",
        "  • figure1_external_roc_ggplot.png / .pdf\n",
        "  • pca_pre_vs_post_combined.png\n",
        "  • biomarker_discovery_outputs.rds\n")

sessionInfo()
