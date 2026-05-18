############################################################
## Generate Figure 2 synthetic datasets:
## sequencing-depth benchmark
##
## Design:
## - ngene = 500 HVGs
## - 4 balanced cell types
## - initial singlets = 2000 total = 500 per cell type
## - target singlet UMI depths = 500, 1000, 2000, 4000, 8000
## - doublet rate = 20%
## - random doublet pairing
## - parent singlets used to form doublets are removed
##
## Key control:
## - The same base scDesign3 singlet profiles are used for all depths.
## - The same parent-pair design is used for all depths.
## - Therefore, only UMI depth changes across datasets.
############################################################


############################################################
## 0. Working directory and packages
############################################################

## Optional: change this to your working directory
setwd("/media/subunit/16T/ICC_scanpy/yyztest/simulate")

library(scDesign3)
library(SingleCellExperiment)
library(DuoClustering2018)
library(scran)
library(scuttle)
library(S4Vectors)
library(dplyr)
library(tibble)
library(Matrix)

set.seed(123)

out_dir <- "simulated_doublet_datasets_Figure2_depth"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


############################################################
## 1. Global parameters
############################################################

ngene <- 500

n_cell_types <- 4
singlets_per_celltype <- 500
initial_singlets_total <- n_cell_types * singlets_per_celltype

target_doublet_rate <- 0.20

target_umi_depths <- c(500, 1000, 2000, 4000, 8000)

cat("Figure 2 sequencing-depth benchmark design:\n")
cat("ngene =", ngene, "\n")
cat("n_cell_types =", n_cell_types, "\n")
cat("singlets_per_celltype =", singlets_per_celltype, "\n")
cat("initial_singlets_total =", initial_singlets_total, "\n")
cat("target_doublet_rate =", target_doublet_rate, "\n")
cat("target_umi_depths =", paste(target_umi_depths, collapse = ", "), "\n")


############################################################
## 2. Read reference data and select 4 cell types
############################################################

sce_ref <- get("sce_filteredExpr10_Zhengmix4eq")(metadata = FALSE)

## True cell-type labels in this dataset are stored as 'phenoid'
colData(sce_ref)$cell_type <- factor(colData(sce_ref)$phenoid)

cat("\nAvailable reference cell types:\n")
print(table(colData(sce_ref)$cell_type))

## Current reference contains exactly 4 cell types:
## b.cells, cd14.monocytes, naive.cytotoxic, regulatory.t
## To keep Figure 2 consistent with Figure 1, use all 4 cell types.
tab_ct <- sort(table(colData(sce_ref)$cell_type), decreasing = TRUE)
selected_ct <- names(tab_ct)[seq_len(min(n_cell_types, length(tab_ct)))]

if (length(selected_ct) != n_cell_types) {
  stop("Could not select exactly 4 cell types. Please check reference data.")
}

cat("\nSelected 4 cell types:\n")
print(selected_ct)

sce_ref <- sce_ref[, colData(sce_ref)$cell_type %in% selected_ct]
colData(sce_ref)$cell_type <- droplevels(factor(colData(sce_ref)$cell_type))

cat("\nReference after selecting 4 cell types:\n")
print(dim(sce_ref))
print(table(colData(sce_ref)$cell_type))


############################################################
## 3. Select top 500 highly variable genes
############################################################

logcounts(sce_ref) <- log1p(counts(sce_ref))

gene_var <- modelGeneVar(sce_ref)
chosen_genes <- getTopHVGs(gene_var, n = min(ngene, nrow(sce_ref)))

sce_ref <- sce_ref[chosen_genes, ]

cat("\nReference after HVG selection:\n")
print(dim(sce_ref))


############################################################
## 4. Fit scDesign3 model once
############################################################

example_data <- construct_data(
  sce = sce_ref,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "cell_type"
)

example_marginal <- fit_marginal(
  data = example_data,
  predictor = "gene",
  mu_formula = "cell_type",
  sigma_formula = "cell_type",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE
)

set.seed(123)
example_copula <- fit_copula(
  sce = sce_ref,
  assay_use = "counts",
  marginal_list = example_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 4,
  input_data = example_data$dat
)


############################################################
## 5. Helper: create balanced 4-cell-type covariate
############################################################

make_balanced_covariate <- function(cell_types, n_per_type) {
  new_cov <- data.frame(
    cell_type = rep(cell_types, each = n_per_type),
    stringsAsFactors = FALSE
  )
  
  new_cov$cell_type <- factor(new_cov$cell_type, levels = cell_types)
  new_cov$corr_group <- new_cov$cell_type
  
  return(new_cov)
}


############################################################
## 6. Helper: resample each cell to fixed UMI depth
############################################################

resample_columns_to_library <- function(mat, target_umi = 2000, seed = 1) {
  set.seed(seed)
  
  mat <- as.matrix(mat)
  
  out <- matrix(
    0L,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  lib <- colSums(mat)
  
  global_prob <- rowSums(mat)
  if (sum(global_prob) <= 0) {
    stop("Global expression is zero; cannot resample columns.")
  }
  global_prob <- global_prob / sum(global_prob)
  
  for (j in seq_len(ncol(mat))) {
    if (lib[j] > 0) {
      prob <- mat[, j] / lib[j]
    } else {
      warning("A zero-library cell was found. Using global gene probabilities for resampling.")
      prob <- global_prob
    }
    
    out[, j] <- as.integer(rmultinom(n = 1, size = target_umi, prob = prob))
  }
  
  return(as(out, "dgCMatrix"))
}


############################################################
## 7. Helper: simulate one base 2000-singlet pool with scDesign3
############################################################

simulate_base_initial_singlets <- function(new_covariate, seed = 1001) {
  set.seed(seed)
  
  example_para <- extract_para(
    sce = sce_ref,
    marginal_list = example_marginal,
    n_cores = 1,
    family_use = "nb",
    new_covariate = new_covariate,
    data = example_data$dat
  )
  
  set.seed(seed)
  new_counts <- simu_new(
    sce = sce_ref,
    mean_mat = example_para$mean_mat,
    sigma_mat = example_para$sigma_mat,
    zero_mat = example_para$zero_mat,
    quantile_mat = NULL,
    copula_list = example_copula$copula_list,
    n_cores = 1,
    family_use = "nb",
    input_data = example_data$dat,
    new_covariate = new_covariate,
    important_feature = example_copula$important_feature,
    filtered_gene = example_data$filtered_gene
  )
  
  new_counts <- as(new_counts, "dgCMatrix")
  
  singlet_meta <- DataFrame(
    cell_type = as.character(new_covariate$cell_type),
    true_status = "singlet",
    doublet_class = "singlet",
    parent1_type = as.character(new_covariate$cell_type),
    parent2_type = NA_character_,
    initial_singlet_id = paste0("initial_singlet_", seq_len(ncol(new_counts)))
  )
  
  sce <- SingleCellExperiment(
    assays = list(counts = new_counts),
    colData = singlet_meta
  )
  
  colnames(sce) <- singlet_meta$initial_singlet_id
  logcounts(sce) <- log1p(counts(sce))
  
  metadata(sce)$ngene <- nrow(sce)
  metadata(sce)$initial_singlets_total <- ncol(sce)
  metadata(sce)$cell_types <- unique(as.character(singlet_meta$cell_type))
  metadata(sce)$note <- "Base scDesign3 singlet pool before UMI-depth resampling"
  
  return(sce)
}


############################################################
## 8. Helper: generate fixed parent-pair design
############################################################

make_parent_pair_design <- function(
    sce_singlet_pool,
    doublet_rate = 0.20,
    seed = 2020
) {
  set.seed(seed)
  
  X <- ncol(sce_singlet_pool)
  Y <- doublet_rate
  
  ## Benchmark-style formula:
  ## D = X * Y / (1 + Y)
  ## final droplets = X - D
  ## observed doublet rate = D / (X - D)
  n_doublets <- round(X * Y / (1 + Y))
  
  if (2 * n_doublets > X) {
    stop("Too many doublets requested: 2 * n_doublets exceeds initial singlet pool.")
  }
  
  all_cells <- seq_len(X)
  
  ## Select 2D parent singlets without replacement.
  parent_idx <- sample(all_cells, size = 2 * n_doublets, replace = FALSE)
  
  ## Randomly shuffle and split into parent1 and parent2.
  parent_idx <- sample(parent_idx, size = length(parent_idx), replace = FALSE)
  
  p1 <- parent_idx[seq_len(n_doublets)]
  p2 <- parent_idx[(n_doublets + 1):(2 * n_doublets)]
  
  used_parent_idx <- c(p1, p2)
  remaining_idx <- setdiff(all_cells, used_parent_idx)
  
  design <- list(
    p1 = p1,
    p2 = p2,
    used_parent_idx = used_parent_idx,
    remaining_idx = remaining_idx,
    n_doublets = n_doublets,
    n_remaining_singlets = length(remaining_idx),
    initial_singlets_total = X,
    target_doublet_rate = doublet_rate,
    observed_doublet_rate = n_doublets / (length(remaining_idx) + n_doublets)
  )
  
  return(design)
}


############################################################
## 9. Helper: create final dataset from one UMI-resampled singlet pool
############################################################

create_depth_dataset_from_design <- function(
    sce_singlet_pool_depth,
    pair_design,
    target_umi_singlet,
    scenario_name
) {
  mat <- counts(sce_singlet_pool_depth)
  cd <- as.data.frame(colData(sce_singlet_pool_depth))
  
  p1 <- pair_design$p1
  p2 <- pair_design$p2
  remaining_idx <- pair_design$remaining_idx
  n_doublets <- pair_design$n_doublets
  
  parent1_type <- as.character(cd$cell_type[p1])
  parent2_type <- as.character(cd$cell_type[p2])
  
  doublet_class <- ifelse(parent1_type == parent2_type, "homotypic", "heterotypic")
  
  ## Synthetic doublet = parent1 count vector + parent2 count vector.
  ## If each parent singlet has target_umi_singlet counts,
  ## each doublet has approximately 2 * target_umi_singlet counts.
  doublet_counts <- mat[, p1, drop = FALSE] + mat[, p2, drop = FALSE]
  colnames(doublet_counts) <- paste0(scenario_name, "_doublet_", seq_len(n_doublets))
  
  singlet_counts <- mat[, remaining_idx, drop = FALSE]
  colnames(singlet_counts) <- paste0(scenario_name, "_singlet_", seq_along(remaining_idx))
  
  singlet_meta <- data.frame(
    cell_type = as.character(cd$cell_type[remaining_idx]),
    true_status = "singlet",
    doublet_class = "singlet",
    parent1_type = as.character(cd$cell_type[remaining_idx]),
    parent2_type = NA_character_,
    initial_singlet_id = cd$initial_singlet_id[remaining_idx],
    parent1_initial_id = NA_character_,
    parent2_initial_id = NA_character_,
    stringsAsFactors = FALSE
  )
  
  doublet_meta <- data.frame(
    cell_type = "doublet",
    true_status = "doublet",
    doublet_class = doublet_class,
    parent1_type = parent1_type,
    parent2_type = parent2_type,
    initial_singlet_id = NA_character_,
    parent1_initial_id = cd$initial_singlet_id[p1],
    parent2_initial_id = cd$initial_singlet_id[p2],
    stringsAsFactors = FALSE
  )
  
  combined_counts <- cbind(singlet_counts, doublet_counts)
  combined_meta <- rbind(singlet_meta, doublet_meta)
  
  sce_out <- SingleCellExperiment(
    assays = list(counts = as(combined_counts, "dgCMatrix")),
    colData = S4Vectors::DataFrame(combined_meta)
  )
  
  colnames(sce_out) <- c(colnames(singlet_counts), colnames(doublet_counts))
  logcounts(sce_out) <- log1p(counts(sce_out))
  
  metadata(sce_out)$scenario <- scenario_name
  metadata(sce_out)$target_umi_singlet <- target_umi_singlet
  metadata(sce_out)$target_doublet_rate <- pair_design$target_doublet_rate
  metadata(sce_out)$observed_doublet_rate <- pair_design$observed_doublet_rate
  metadata(sce_out)$initial_singlets_total <- pair_design$initial_singlets_total
  metadata(sce_out)$n_doublets_generated <- pair_design$n_doublets
  metadata(sce_out)$n_parent_singlets_removed <- 2 * pair_design$n_doublets
  metadata(sce_out)$n_remaining_singlets <- pair_design$n_remaining_singlets
  metadata(sce_out)$n_cell_types <- length(unique(cd$cell_type))
  metadata(sce_out)$doublet_generation <- "random_pairing_parent_singlets_removed"
  metadata(sce_out)$doublet_count_rule <- "parent1_counts_plus_parent2_counts"
  metadata(sce_out)$control_note <- "Same base singlet pool and same parent-pair design used across all UMI depths"
  
  return(sce_out)
}


############################################################
## 10. Helper: sanity check
############################################################

check_simulated_dataset <- function(sce) {
  mat <- counts(sce)
  cd <- as.data.frame(colData(sce))
  
  singlet_idx <- cd$true_status == "singlet"
  doublet_idx <- cd$true_status == "doublet"
  
  data.frame(
    scenario = metadata(sce)$scenario,
    n_genes = nrow(sce),
    n_cells = ncol(sce),
    n_cell_types = metadata(sce)$n_cell_types,
    initial_singlets_total = metadata(sce)$initial_singlets_total,
    n_singlets = sum(singlet_idx),
    n_doublets = sum(doublet_idx),
    target_doublet_rate = metadata(sce)$target_doublet_rate,
    observed_doublet_rate = mean(doublet_idx),
    target_umi_singlet = metadata(sce)$target_umi_singlet,
    n_parent_singlets_removed = metadata(sce)$n_parent_singlets_removed,
    object_size_MB = as.numeric(object.size(sce)) / 1024^2,
    zero_fraction = 1 - Matrix::nnzero(mat) / (nrow(mat) * ncol(mat)),
    median_UMI_singlet = median(Matrix::colSums(mat)[singlet_idx]),
    median_UMI_doublet = median(Matrix::colSums(mat)[doublet_idx]),
    median_genes_singlet = median(Matrix::colSums(mat > 0)[singlet_idx]),
    median_genes_doublet = median(Matrix::colSums(mat > 0)[doublet_idx]),
    n_homotypic_doublets = sum(cd$doublet_class == "homotypic", na.rm = TRUE),
    n_heterotypic_doublets = sum(cd$doublet_class == "heterotypic", na.rm = TRUE),
    heterotypic_doublet_fraction = sum(cd$doublet_class == "heterotypic", na.rm = TRUE) /
      sum(doublet_idx),
    stringsAsFactors = FALSE
  )
}


############################################################
## 11. Generate base initial singlet pool
############################################################

new_cov <- make_balanced_covariate(
  cell_types = levels(colData(sce_ref)$cell_type),
  n_per_type = singlets_per_celltype
)

cat("\nNew covariate table for synthetic singlets:\n")
print(table(new_cov$cell_type))

base_singlet_sce <- simulate_base_initial_singlets(
  new_covariate = new_cov,
  seed = 1001
)

cat("\nBase synthetic singlet pool before UMI-depth resampling:\n")
print(dim(base_singlet_sce))
print(table(colData(base_singlet_sce)$cell_type))
cat("Median raw base singlet UMI before resampling:\n")
print(median(Matrix::colSums(counts(base_singlet_sce))))

saveRDS(
  base_singlet_sce,
  file = file.path(out_dir, "Depth_base_scDesign3_2000_singlets_4celltypes_unfixedUMI.rds")
)


############################################################
## 12. Generate one fixed parent-pair design
############################################################

pair_design <- make_parent_pair_design(
  sce_singlet_pool = base_singlet_sce,
  doublet_rate = target_doublet_rate,
  seed = 2020
)

cat("\nFixed parent-pair design:\n")
print(pair_design[c(
  "initial_singlets_total",
  "n_doublets",
  "n_remaining_singlets",
  "target_doublet_rate",
  "observed_doublet_rate"
)])

saveRDS(
  pair_design,
  file = file.path(out_dir, "Depth_fixed_parent_pair_design_D20.rds")
)


############################################################
## 13. Generate five sequencing-depth datasets
############################################################

scenario_tbl <- tibble::tibble(
  scenario = c(
    "UMI0.5k_D20",
    "UMI1k_D20",
    "UMI2k_D20",
    "UMI4k_D20",
    "UMI8k_D20"
  ),
  target_umi_singlet = target_umi_depths,
  seed = c(500, 1000, 2000, 4000, 8000)
)

cat("\nScenario table:\n")
print(scenario_tbl)

simulated_datasets <- list()
sanity_list <- list()

for (i in seq_len(nrow(scenario_tbl))) {
  sc <- scenario_tbl[i, ]
  
  cat("\nGenerating ", sc$scenario, " ...\n", sep = "")
  
  ## Resample the same base singlet pool to the target UMI depth.
  depth_counts <- resample_columns_to_library(
    mat = counts(base_singlet_sce),
    target_umi = sc$target_umi_singlet,
    seed = sc$seed + 10
  )
  
  depth_singlet_sce <- base_singlet_sce
  assays(depth_singlet_sce)$counts <- depth_counts
  logcounts(depth_singlet_sce) <- log1p(counts(depth_singlet_sce))
  metadata(depth_singlet_sce)$target_umi_singlet <- sc$target_umi_singlet
  
  ## Create final contaminated dataset using the same parent-pair design.
  final_sce <- create_depth_dataset_from_design(
    sce_singlet_pool_depth = depth_singlet_sce,
    pair_design = pair_design,
    target_umi_singlet = sc$target_umi_singlet,
    scenario_name = sc$scenario
  )
  
  check_row <- check_simulated_dataset(final_sce)
  print(check_row)
  
  if (abs(check_row$observed_doublet_rate - check_row$target_doublet_rate) > 0.005) {
    warning("Observed doublet rate differs from target rate by >0.5% in ", sc$scenario)
  }
  
  if (check_row$median_UMI_singlet != sc$target_umi_singlet) {
    warning("Median singlet UMI is not exactly target UMI in ", sc$scenario)
  }
  
  if (check_row$median_UMI_doublet <= check_row$median_UMI_singlet) {
    warning("Doublets do not have higher median UMI than singlets in ", sc$scenario)
  }
  
  if (check_row$median_genes_doublet <= check_row$median_genes_singlet) {
    warning("Doublets do not have more detected genes than singlets in ", sc$scenario)
  }
  
  simulated_datasets[[sc$scenario]] <- final_sce
  sanity_list[[sc$scenario]] <- check_row
  
  saveRDS(
    final_sce,
    file = file.path(out_dir, paste0(sc$scenario, ".rds"))
  )
}


############################################################
## 14. Save summary tables
############################################################

dataset_sanity_check <- do.call(rbind, sanity_list)

write.csv(
  dataset_sanity_check,
  file = file.path(out_dir, "Figure2_depth_dataset_sanity_check.csv"),
  row.names = FALSE
)

dataset_summary <- dataset_sanity_check %>%
  dplyr::select(
    scenario,
    n_genes,
    n_cell_types,
    initial_singlets_total,
    n_singlets,
    n_doublets,
    n_cells,
    target_doublet_rate,
    observed_doublet_rate,
    target_umi_singlet,
    n_parent_singlets_removed,
    median_UMI_singlet,
    median_UMI_doublet,
    median_genes_singlet,
    median_genes_doublet,
    n_homotypic_doublets,
    n_heterotypic_doublets,
    heterotypic_doublet_fraction,
    zero_fraction
  )

write.csv(
  dataset_summary,
  file = file.path(out_dir, "Figure2_depth_dataset_summary.csv"),
  row.names = FALSE
)

cat("\nAll Figure 2 sequencing-depth datasets saved to: ", out_dir, "\n", sep = "")
cat("\nDataset summary:\n")
print(dataset_summary)

sessionInfo()