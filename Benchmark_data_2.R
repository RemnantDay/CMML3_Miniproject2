############################################################
## Figure 2 sequencing-depth benchmark
##
## Input datasets:
## simulated_doublet_datasets_Figure2_depth/
##   UMI0.5k_D20.rds
##   UMI1k_D20.rds
##   UMI2k_D20.rds
##   UMI4k_D20.rds
##   UMI8k_D20.rds
##
## Methods:
## - scDblFinder_random
## - scds_cxds
## - scds_bcds
## - scds_hybrid
## - DoubletFinder
##
## Metrics:
## - AUPRC
## - AUROC
## - runtime_sec
## - homotypic_recall
## - heterotypic_recall
## - matched_recovery
## - FDP
############################################################


############################################################
## 0. Packages
############################################################
setwd("/media/subunit/16T/ICC_scanpy/yyztest/simulate")
library(SingleCellExperiment)
library(Matrix)

library(scDblFinder)
library(scds)

library(Seurat)
library(DoubletFinder)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)

library(pROC)
library(PRROC)

set.seed(123)


############################################################
## 1. Paths and scenario table
############################################################

data_dir <- "simulated_doublet_datasets_Figure2_depth"
out_dir <- file.path(data_dir, "benchmark_Figure2_depth")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scenario_tbl <- tibble::tribble(
  ~scenario,        ~target_umi_singlet,
  "UMI0.5k_D20",     500,
  "UMI1k_D20",      1000,
  "UMI2k_D20",      2000,
  "UMI4k_D20",      4000,
  "UMI8k_D20",      8000
) %>%
  mutate(
    file = file.path(data_dir, paste0(scenario, ".rds")),
    target_umi_k = target_umi_singlet / 1000
  )

print(scenario_tbl)

if (!all(file.exists(scenario_tbl$file))) {
  stop("Some Figure 2 RDS files are missing. Please check data_dir and filenames.")
}


############################################################
## 2. Helper functions
############################################################

clean_score <- function(x) {
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    stop("All scores are NA.")
  }
  
  finite_idx <- is.finite(x)
  if (!any(finite_idx)) {
    stop("No finite scores available.")
  }
  
  min_finite <- min(x[finite_idx], na.rm = TRUE)
  max_finite <- max(x[finite_idx], na.rm = TRUE)
  
  x[is.na(x)] <- min_finite - 1e-8
  x[x == Inf] <- max_finite + 1e-8
  x[x == -Inf] <- min_finite - 1e-8
  
  return(x)
}


make_rate_matched_calls <- function(score, true_status) {
  score <- clean_score(score)
  
  n_true_doublets <- sum(true_status == "doublet")
  
  ord <- order(score, decreasing = TRUE)
  
  pred <- rep("singlet", length(score))
  pred[ord[seq_len(n_true_doublets)]] <- "doublet"
  
  return(pred)
}


calc_metrics <- function(score, true_status, doublet_class, pred_class) {
  score <- clean_score(score)
  y_true <- ifelse(true_status == "doublet", 1, 0)
  
  ## AUROC: score higher = more doublet-like
  auroc <- as.numeric(
    pROC::roc(
      response = y_true,
      predictor = score,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )$auc
  )
  
  ## AUPRC: positive class = true doublets
  auprc <- PRROC::pr.curve(
    scores.class0 = score[y_true == 1],
    scores.class1 = score[y_true == 0],
    curve = FALSE
  )$auc.integral
  
  tp <- sum(pred_class == "doublet" & true_status == "doublet")
  fp <- sum(pred_class == "doublet" & true_status == "singlet")
  tn <- sum(pred_class == "singlet" & true_status == "singlet")
  fn <- sum(pred_class == "singlet" & true_status == "doublet")
  
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  f1 <- ifelse(
    is.na(precision) | is.na(recall) | (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  fdp <- ifelse((tp + fp) == 0, NA_real_, fp / (tp + fp))
  tnr <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  fpr <- ifelse((tn + fp) == 0, NA_real_, fp / (tn + fp))
  accuracy <- (tp + tn) / (tp + fp + tn + fn)
  
  homotypic_idx <- true_status == "doublet" & doublet_class == "homotypic"
  heterotypic_idx <- true_status == "doublet" & doublet_class == "heterotypic"
  
  homotypic_recall <- ifelse(
    sum(homotypic_idx) == 0,
    NA_real_,
    sum(pred_class == "doublet" & homotypic_idx) / sum(homotypic_idx)
  )
  
  heterotypic_recall <- ifelse(
    sum(heterotypic_idx) == 0,
    NA_real_,
    sum(pred_class == "doublet" & heterotypic_idx) / sum(heterotypic_idx)
  )
  
  tibble(
    AUROC = auroc,
    AUPRC = auprc,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    FDP = fdp,
    TNR = tnr,
    FPR = fpr,
    TP = tp,
    FP = fp,
    TN = tn,
    FN = fn,
    matched_recovery = recall,
    homotypic_recall = homotypic_recall,
    heterotypic_recall = heterotypic_recall,
    n_true_doublets = sum(true_status == "doublet"),
    n_pred_doublets = sum(pred_class == "doublet"),
    n_true_singlets = sum(true_status == "singlet")
  )
}


make_call_outcome <- function(true_status, doublet_class, pred_class) {
  dplyr::case_when(
    pred_class != "doublet" ~ NA_character_,
    true_status == "singlet" ~ "False-positive singlet",
    true_status == "doublet" & doublet_class == "heterotypic" ~ "True heterotypic doublet",
    true_status == "doublet" & doublet_class == "homotypic" ~ "True homotypic doublet",
    TRUE ~ "Other"
  )
}


evaluate_one_method_on_one_dataset <- function(
    sce,
    score,
    method,
    scenario,
    target_umi_singlet,
    target_umi_k,
    runtime_sec
) {
  meta <- as.data.frame(colData(sce))
  
  true_status <- as.character(meta$true_status)
  doublet_class <- as.character(meta$doublet_class)
  
  score <- clean_score(score)
  
  ## Rate-matched threshold:
  ## predicted doublet number = true doublet number.
  pred_rate_matched <- make_rate_matched_calls(
    score = score,
    true_status = true_status
  )
  
  metrics_tbl <- calc_metrics(
    score = score,
    true_status = true_status,
    doublet_class = doublet_class,
    pred_class = pred_rate_matched
  ) %>%
    mutate(
      scenario = scenario,
      target_umi_singlet = target_umi_singlet,
      target_umi_k = target_umi_k,
      method = method,
      threshold_type = "rate_matched",
      runtime_sec = runtime_sec,
      target_doublet_rate = metadata(sce)$target_doublet_rate,
      observed_doublet_rate = metadata(sce)$observed_doublet_rate,
      .before = 1
    )
  
  pred_tbl <- tibble(
    scenario = scenario,
    target_umi_singlet = target_umi_singlet,
    target_umi_k = target_umi_k,
    method = method,
    cell_id = colnames(sce),
    score = score,
    true_status = true_status,
    doublet_class = doublet_class,
    pred_rate_matched = pred_rate_matched
  ) %>%
    mutate(
      call_outcome = make_call_outcome(
        true_status = true_status,
        doublet_class = doublet_class,
        pred_class = pred_rate_matched
      )
    )
  
  return(list(
    metrics = metrics_tbl,
    predictions = pred_tbl
  ))
}


############################################################
## 3. Method 1: scDblFinder_random
############################################################

metrics_scDblFinder <- list()
pred_scDblFinder <- list()

for (i in seq_len(nrow(scenario_tbl))) {
  sc <- scenario_tbl[i, ]
  
  cat("\n[scDblFinder_random] Running ", sc$scenario, "...\n", sep = "")
  
  sce <- readRDS(sc$file)
  
  runtime <- system.time({
    set.seed(123)
    sce_out <- scDblFinder::scDblFinder(
      sce,
      clusters = FALSE
    )
  })["elapsed"]
  
  score <- colData(sce_out)$scDblFinder.score
  
  res <- evaluate_one_method_on_one_dataset(
    sce = sce,
    score = score,
    method = "scDblFinder_random",
    scenario = sc$scenario,
    target_umi_singlet = sc$target_umi_singlet,
    target_umi_k = sc$target_umi_k,
    runtime_sec = as.numeric(runtime)
  )
  
  metrics_scDblFinder[[sc$scenario]] <- res$metrics
  pred_scDblFinder[[sc$scenario]] <- res$predictions
  
  print(
    res$metrics %>%
      select(
        scenario,
        method,
        target_umi_singlet,
        AUROC,
        AUPRC,
        matched_recovery,
        FDP,
        homotypic_recall,
        heterotypic_recall,
        runtime_sec
      )
  )
}

metrics_scDblFinder_df <- bind_rows(metrics_scDblFinder)
pred_scDblFinder_df <- bind_rows(pred_scDblFinder)

write.csv(
  metrics_scDblFinder_df,
  file = file.path(out_dir, "metrics_scDblFinder_random.csv"),
  row.names = FALSE
)

write.csv(
  pred_scDblFinder_df,
  file = file.path(out_dir, "predictions_scDblFinder_random.csv"),
  row.names = FALSE
)


############################################################
## 4. Methods 2-4: scds_cxds, scds_bcds, scds_hybrid
############################################################

metrics_scds <- list()
pred_scds <- list()

for (i in seq_len(nrow(scenario_tbl))) {
  sc <- scenario_tbl[i, ]
  
  cat("\n[scds] Running ", sc$scenario, "...\n", sep = "")
  
  sce <- readRDS(sc$file)
  
  runtime <- system.time({
    set.seed(123)
    sce_out <- scds::cxds_bcds_hybrid(sce)
  })["elapsed"]
  
  cd_names <- colnames(as.data.frame(colData(sce_out)))
  
  required_cols <- c("cxds_score", "bcds_score", "hybrid_score")
  if (!all(required_cols %in% cd_names)) {
    cat("Available colData columns:\n")
    print(cd_names)
    stop("Expected scds score columns were not found.")
  }
  
  score_list <- list(
    scds_cxds = colData(sce_out)$cxds_score,
    scds_bcds = colData(sce_out)$bcds_score,
    scds_hybrid = colData(sce_out)$hybrid_score
  )
  
  for (m in names(score_list)) {
    res <- evaluate_one_method_on_one_dataset(
      sce = sce,
      score = score_list[[m]],
      method = m,
      scenario = sc$scenario,
      target_umi_singlet = sc$target_umi_singlet,
      target_umi_k = sc$target_umi_k,
      runtime_sec = as.numeric(runtime)
    )
    
    metrics_scds[[paste(sc$scenario, m, sep = "_")]] <- res$metrics
    pred_scds[[paste(sc$scenario, m, sep = "_")]] <- res$predictions
    
    print(
      res$metrics %>%
        select(
          scenario,
          method,
          target_umi_singlet,
          AUROC,
          AUPRC,
          matched_recovery,
          FDP,
          homotypic_recall,
          heterotypic_recall,
          runtime_sec
        )
    )
  }
}

metrics_scds_df <- bind_rows(metrics_scds)
pred_scds_df <- bind_rows(pred_scds)

write.csv(
  metrics_scds_df,
  file = file.path(out_dir, "metrics_scds_cxds_bcds_hybrid.csv"),
  row.names = FALSE
)

write.csv(
  pred_scds_df,
  file = file.path(out_dir, "predictions_scds_cxds_bcds_hybrid.csv"),
  row.names = FALSE
)


############################################################
## 5. Method 5: DoubletFinder wrapper
############################################################

run_DoubletFinder_wrapper <- function(
    sce,
    n_pcs = 20,
    pN = 0.25,
    seed = 123
) {
  set.seed(seed)
  
  mat <- counts(sce)
  true_status <- as.character(colData(sce)$true_status)
  nExp <- sum(true_status == "doublet")
  
  seu <- Seurat::CreateSeuratObject(counts = mat)
  
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  
  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = min(2000, nrow(seu)),
    verbose = FALSE
  )
  
  seu <- Seurat::ScaleData(
    seu,
    features = rownames(seu),
    verbose = FALSE
  )
  
  npcs_use <- min(n_pcs, nrow(seu) - 1, ncol(seu) - 1)
  
  seu <- Seurat::RunPCA(
    seu,
    npcs = npcs_use,
    verbose = FALSE
  )
  
  df_exports <- getNamespaceExports("DoubletFinder")
  
  if ("paramSweep" %in% df_exports) {
    paramSweep_fun <- DoubletFinder::paramSweep
  } else if ("paramSweep_v3" %in% df_exports) {
    paramSweep_fun <- DoubletFinder::paramSweep_v3
  } else {
    stop("Neither paramSweep nor paramSweep_v3 found in DoubletFinder.")
  }
  
  if ("doubletFinder" %in% df_exports) {
    doubletFinder_fun <- DoubletFinder::doubletFinder
  } else if ("doubletFinder_v3" %in% df_exports) {
    doubletFinder_fun <- DoubletFinder::doubletFinder_v3
  } else {
    stop("Neither doubletFinder nor doubletFinder_v3 found in DoubletFinder.")
  }
  
  sweep.res <- paramSweep_fun(
    seu,
    PCs = 1:npcs_use,
    sct = FALSE
  )
  
  sweep.stats <- DoubletFinder::summarizeSweep(
    sweep.res,
    GT = FALSE
  )
  
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  bcmvn <- as.data.frame(bcmvn)
  
  bcmvn$pK_num <- suppressWarnings(as.numeric(as.character(bcmvn$pK)))
  bcmvn$BCmetric_num <- suppressWarnings(as.numeric(as.character(bcmvn$BCmetric)))
  
  bcmvn_valid <- bcmvn[
    is.finite(bcmvn$pK_num) &
      is.finite(bcmvn$BCmetric_num),
  ]
  
  if (nrow(bcmvn_valid) == 0) {
    warning("No valid pK found. Using pK = 0.09 as fallback.")
    best_pK <- 0.09
  } else {
    best_pK <- bcmvn_valid$pK_num[which.max(bcmvn_valid$BCmetric_num)]
  }
  
  best_pK <- as.numeric(best_pK)[1]
  
  message("DoubletFinder selected pK = ", best_pK, "; nExp = ", nExp)
  
  ## Important: reuse.pANN = NULL, not FALSE.
  seu <- doubletFinder_fun(
    seu,
    PCs = 1:npcs_use,
    pN = pN,
    pK = best_pK,
    nExp = nExp,
    reuse.pANN = NULL,
    sct = FALSE
  )
  
  meta_df <- seu@meta.data
  
  pANN_col <- grep("^pANN", colnames(meta_df), value = TRUE)
  class_col <- grep("^DF.classifications", colnames(meta_df), value = TRUE)
  
  if (length(pANN_col) == 0) {
    stop(
      "DoubletFinder pANN column not found. Existing metadata columns are:\n",
      paste(colnames(meta_df), collapse = ", ")
    )
  }
  
  pANN_col <- tail(pANN_col, 1)
  score <- as.numeric(meta_df[[pANN_col]])
  
  if (length(class_col) == 0) {
    default_class <- NA_character_
  } else {
    class_col <- tail(class_col, 1)
    default_class <- ifelse(meta_df[[class_col]] == "Doublet", "doublet", "singlet")
  }
  
  list(
    score = clean_score(score),
    default_class = default_class,
    selected_pK = best_pK,
    nExp = nExp,
    npcs = npcs_use
  )
}


############################################################
## 6. Method 5: DoubletFinder benchmark
############################################################

metrics_DoubletFinder <- list()
pred_DoubletFinder <- list()

for (i in seq_len(nrow(scenario_tbl))) {
  sc <- scenario_tbl[i, ]
  
  cat("\n[DoubletFinder] Running ", sc$scenario, "...\n", sep = "")
  
  sce <- readRDS(sc$file)
  
  runtime <- system.time({
    set.seed(123)
    df_out <- run_DoubletFinder_wrapper(
      sce = sce,
      n_pcs = 20,
      pN = 0.25,
      seed = 123
    )
  })["elapsed"]
  
  score <- df_out$score
  
  res <- evaluate_one_method_on_one_dataset(
    sce = sce,
    score = score,
    method = "DoubletFinder",
    scenario = sc$scenario,
    target_umi_singlet = sc$target_umi_singlet,
    target_umi_k = sc$target_umi_k,
    runtime_sec = as.numeric(runtime)
  )
  
  res$metrics <- res$metrics %>%
    mutate(
      selected_pK = df_out$selected_pK,
      nExp_DoubletFinder = df_out$nExp,
      npcs_DoubletFinder = df_out$npcs
    )
  
  metrics_DoubletFinder[[sc$scenario]] <- res$metrics
  pred_DoubletFinder[[sc$scenario]] <- res$predictions
  
  print(
    res$metrics %>%
      select(
        scenario,
        method,
        target_umi_singlet,
        AUROC,
        AUPRC,
        matched_recovery,
        FDP,
        homotypic_recall,
        heterotypic_recall,
        selected_pK,
        runtime_sec
      )
  )
}

metrics_DoubletFinder_df <- bind_rows(metrics_DoubletFinder)
pred_DoubletFinder_df <- bind_rows(pred_DoubletFinder)

write.csv(
  metrics_DoubletFinder_df,
  file = file.path(out_dir, "metrics_DoubletFinder.csv"),
  row.names = FALSE
)

write.csv(
  pred_DoubletFinder_df,
  file = file.path(out_dir, "predictions_DoubletFinder.csv"),
  row.names = FALSE
)


############################################################
## 7. Combine all benchmark outputs
############################################################

metric_files <- c(
  file.path(out_dir, "metrics_scDblFinder_random.csv"),
  file.path(out_dir, "metrics_scds_cxds_bcds_hybrid.csv"),
  file.path(out_dir, "metrics_DoubletFinder.csv")
)

prediction_files <- c(
  file.path(out_dir, "predictions_scDblFinder_random.csv"),
  file.path(out_dir, "predictions_scds_cxds_bcds_hybrid.csv"),
  file.path(out_dir, "predictions_DoubletFinder.csv")
)

if (!all(file.exists(metric_files))) {
  stop("Some metric files are missing. Run all method sections first.")
}

if (!all(file.exists(prediction_files))) {
  stop("Some prediction files are missing. Run all method sections first.")
}

metrics_all <- metric_files %>%
  lapply(read.csv) %>%
  bind_rows()

predictions_all <- prediction_files %>%
  lapply(read.csv) %>%
  bind_rows()

method_order <- c(
  "scDblFinder_random",
  "scds_cxds",
  "scds_bcds",
  "scds_hybrid",
  "DoubletFinder"
)

scenario_order <- scenario_tbl$scenario

metrics_all <- metrics_all %>%
  mutate(
    method = factor(method, levels = method_order),
    scenario = factor(scenario, levels = scenario_order),
    target_umi_singlet = as.numeric(target_umi_singlet),
    target_umi_k = as.numeric(target_umi_k)
  ) %>%
  arrange(target_umi_singlet, method)

predictions_all <- predictions_all %>%
  mutate(
    method = factor(method, levels = method_order),
    scenario = factor(scenario, levels = scenario_order),
    target_umi_singlet = as.numeric(target_umi_singlet),
    target_umi_k = as.numeric(target_umi_k)
  )

write.csv(
  metrics_all,
  file = file.path(out_dir, "Figure2_depth_all_metrics.csv"),
  row.names = FALSE
)

write.csv(
  predictions_all,
  file = file.path(out_dir, "Figure2_depth_all_cell_level_predictions.csv"),
  row.names = FALSE
)

cat("\nCombined Figure 2 depth benchmark metrics:\n")
print(
  metrics_all %>%
    select(
      scenario,
      target_umi_singlet,
      method,
      AUROC,
      AUPRC,
      matched_recovery,
      FDP,
      homotypic_recall,
      heterotypic_recall,
      runtime_sec
    )
)


############################################################
## 8. Prepare plot-ready tables
############################################################

## AUPRC and AUROC
fig2_score_df <- metrics_all %>%
  select(
    scenario,
    target_umi_singlet,
    target_umi_k,
    method,
    AUPRC,
    AUROC
  ) %>%
  pivot_longer(
    cols = c(AUPRC, AUROC),
    names_to = "metric",
    values_to = "value"
  )

write.csv(
  fig2_score_df,
  file = file.path(out_dir, "Figure2_depth_AUPRC_AUROC_plotdata.csv"),
  row.names = FALSE
)

## Runtime
fig2_runtime_df <- metrics_all %>%
  select(
    scenario,
    target_umi_singlet,
    target_umi_k,
    method,
    runtime_sec
  )

write.csv(
  fig2_runtime_df,
  file = file.path(out_dir, "Figure2_depth_runtime_plotdata.csv"),
  row.names = FALSE
)

## Homotypic / heterotypic recall
fig2_class_recall_df <- metrics_all %>%
  select(
    scenario,
    target_umi_singlet,
    target_umi_k,
    method,
    homotypic_recall,
    heterotypic_recall
  ) %>%
  pivot_longer(
    cols = c(homotypic_recall, heterotypic_recall),
    names_to = "doublet_class",
    values_to = "recall"
  ) %>%
  mutate(
    doublet_class = recode(
      doublet_class,
      homotypic_recall = "Homotypic doublets",
      heterotypic_recall = "Heterotypic doublets"
    )
  )

write.csv(
  fig2_class_recall_df,
  file = file.path(out_dir, "Figure2_depth_homotypic_heterotypic_recall_plotdata.csv"),
  row.names = FALSE
)

## Summary table for report
compact_metrics <- metrics_all %>%
  select(
    scenario,
    target_umi_singlet,
    target_umi_k,
    method,
    AUROC,
    AUPRC,
    matched_recovery,
    FDP,
    homotypic_recall,
    heterotypic_recall,
    runtime_sec,
    TP,
    FP,
    TN,
    FN
  )

write.csv(
  compact_metrics,
  file = file.path(out_dir, "Figure2_depth_compact_metrics_table.csv"),
  row.names = FALSE
)

cat("\nAll Figure 2 benchmark outputs saved to:\n")
cat(out_dir, "\n")

sessionInfo()