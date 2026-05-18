############################################################
## Figure 1 benchmark: common setup
############################################################

setwd("/media/subunit/16T/ICC_scanpy/yyztest/simulate")

library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)

library(scDblFinder)
library(scds)
library(scran)

library(pROC)
library(PRROC)

set.seed(123)

data_dir <- "simulated_doublet_datasets_Figure1_4CT"
out_dir <- file.path(data_dir, "benchmark_Figure1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scenario_tbl <- tibble::tribble(
  ~scenario,              ~target_doublet_rate,
  "B5_randomDoublets",     0.05,
  "B10_randomDoublets",    0.10,
  "B20_randomDoublets",    0.20,
  "B30_randomDoublets",    0.30,
  "B40_randomDoublets",    0.40
) %>%
  mutate(
    file = file.path(data_dir, paste0(scenario, ".rds")),
    target_doublet_rate_pct = target_doublet_rate * 100
  )

print(scenario_tbl)

if (!all(file.exists(scenario_tbl$file))) {
  stop("Some Figure 1 RDS files are missing. Check data_dir and filenames.")
}

############################################################
## Helper functions
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
  recall    <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  f1        <- ifelse(
    is.na(precision) | is.na(recall) | (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  fdp <- ifelse((tp + fp) == 0, NA_real_, fp / (tp + fp))
  tnr <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  fpr <- ifelse((tn + fp) == 0, NA_real_, fp / (tn + fp))
  
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
    target_doublet_rate,
    target_doublet_rate_pct,
    runtime_sec
) {
  meta <- as.data.frame(colData(sce))
  
  true_status <- as.character(meta$true_status)
  doublet_class <- as.character(meta$doublet_class)
  
  score <- clean_score(score)
  
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
      target_doublet_rate = target_doublet_rate,
      target_doublet_rate_pct = target_doublet_rate_pct,
      method = method,
      runtime_sec = runtime_sec,
      .before = 1
    )
  
  pred_tbl <- tibble(
    scenario = scenario,
    target_doublet_rate = target_doublet_rate,
    target_doublet_rate_pct = target_doublet_rate_pct,
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
## Method 1: scDblFinder_random
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
    target_doublet_rate = sc$target_doublet_rate,
    target_doublet_rate_pct = sc$target_doublet_rate_pct,
    runtime_sec = as.numeric(runtime)
  )
  
  metrics_scDblFinder[[sc$scenario]] <- res$metrics
  pred_scDblFinder[[sc$scenario]] <- res$predictions
  
  print(res$metrics %>% select(scenario, method, AUROC, AUPRC, Recall, FDP, homotypic_recall, heterotypic_recall))
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

## Method 1 successfully finished.

############################################################
## Methods 2-4: scds_cxds, scds_bcds, scds_hybrid
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
  print(cd_names)
  
  required_cols <- c("cxds_score", "bcds_score", "hybrid_score")
  if (!all(required_cols %in% cd_names)) {
    stop("Expected scds score columns not found. Check printed colData column names.")
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
      target_doublet_rate = sc$target_doublet_rate,
      target_doublet_rate_pct = sc$target_doublet_rate_pct,
      runtime_sec = as.numeric(runtime)
    )
    
    metrics_scds[[paste(sc$scenario, m, sep = "_")]] <- res$metrics
    pred_scds[[paste(sc$scenario, m, sep = "_")]] <- res$predictions
    
    print(res$metrics %>% select(scenario, method, AUROC, AUPRC, Recall, FDP, homotypic_recall, heterotypic_recall))
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
## Method 5: DoubletFinder
############################################################

metrics_DoubletFinder <- list()
pred_DoubletFinder <- list()
library(Seurat)
library(DoubletFinder)
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
  
  ## Convert to Seurat object
  seu <- Seurat::CreateSeuratObject(counts = mat)
  
  ## Standard preprocessing
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
  
  ## Choose correct function names depending on DoubletFinder version
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
  
  ## Parameter sweep to choose pK
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
  
  ## Run DoubletFinder
  ## IMPORTANT:
  ## use reuse.pANN = NULL, not FALSE.
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
  
  ## If multiple pANN columns exist, use the newest one.
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
    target_doublet_rate = sc$target_doublet_rate,
    target_doublet_rate_pct = sc$target_doublet_rate_pct,
    runtime_sec = as.numeric(runtime)
  )
  
  ## Add selected pK and nExp to metrics table
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
        AUROC,
        AUPRC,
        Recall,
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
