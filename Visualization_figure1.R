setwd("/media/subunit/16T/ICC_scanpy/yyztest/simulate")
############################################################
## Figure 1 visualization
## 6 plots in a 2 x 3 layout
##
## Inputs:
## - metrics_scDblFinder_random.csv
## - metrics_scds_cxds_bcds_hybrid.csv
## - metrics_DoubletFinder.csv
##
## Outputs:
## - 5 line plots
## - 1 outperforming-count heatmap
## - combined 2x3 figure
############################################################

############################
## 0. Packages
############################
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(readr)
library(grid)

############################
## 1. Paths
############################
data_dir <- "simulated_doublet_datasets_Figure1_4CT"
benchmark_dir <- file.path(data_dir, "benchmark_Figure1")
plot_dir <- file.path(benchmark_dir, "Figure1_visualization_updated")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

metric_files <- c(
  file.path(benchmark_dir, "metrics_scDblFinder_random.csv"),
  file.path(benchmark_dir, "metrics_scds_cxds_bcds_hybrid.csv"),
  file.path(benchmark_dir, "metrics_DoubletFinder.csv")
)

if (!all(file.exists(metric_files))) {
  stop("Some metrics CSV files are missing. Please check benchmark_dir.")
}

############################
## 2. Read and combine data
############################
metrics_all <- metric_files %>%
  lapply(read.csv) %>%
  bind_rows()

cat("Columns in metrics_all:\n")
print(colnames(metrics_all))

required_cols <- c(
  "scenario",
  "method",
  "target_doublet_rate_pct",
  "AUPRC",
  "AUROC",
  "homotypic_recall",
  "heterotypic_recall",
  "runtime_sec"
)

missing_cols <- setdiff(required_cols, colnames(metrics_all))
if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing:\n",
    paste(missing_cols, collapse = ", ")
  )
}

############################
## 3. Factor levels and colors
############################
method_order <- c(
  "scDblFinder_random",
  "scds_cxds",
  "scds_bcds",
  "scds_hybrid",
  "DoubletFinder"
)

scenario_order <- c(
  "B5_randomDoublets",
  "B10_randomDoublets",
  "B20_randomDoublets",
  "B30_randomDoublets",
  "B40_randomDoublets"
)

method_colors <- c(
  "scDblFinder_random" = "#E31A1B",
  "scds_cxds"          = "#367db7",
  "scds_bcds"          = "#4FAD4A",
  "scds_hybrid"        = "#974DA1",
  "DoubletFinder"      = "#FD7F00"
)

method_labels <- c(
  "scDblFinder_random" = "scDblFinder",
  "scds_cxds"          = "cxds",
  "scds_bcds"          = "bcds",
  "scds_hybrid"        = "hybrid",
  "DoubletFinder"      = "DoubletFinder"
)

metrics_all <- metrics_all %>%
  mutate(
    method = factor(method, levels = method_order),
    scenario = factor(scenario, levels = scenario_order),
    target_doublet_rate_pct = as.numeric(target_doublet_rate_pct),
    run_time = runtime_sec
  ) %>%
  arrange(target_doublet_rate_pct, method)

############################
## 4. Theme settings
############################
base_theme <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size=9),
    plot.subtitle = element_text(hjust = 0.5,size=9),
    axis.title = element_text(face = "bold",size=8),
    axis.text = element_text(color = "black",size=8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.15),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.15),
    legend.title = element_text(face = "bold",size=8),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold",size=9)
  )

## Legend inside plot, on the right side
legend_inside_right <- theme(
  legend.position = 'right',
  legend.justification = 'center',
  legend.background = element_blank(),
  legend.key.height = unit(0.35, "cm"),
  legend.key.width = unit(0.45, "cm")
)

############################
## 5. Helper function for line plots
############################
make_line_plot <- function(df, ycol, ylab, title_text,
                           show_legend = FALSE,
                           y_limits = NULL,
                           y_breaks = waiver()) {
  
  p <- ggplot(
    df,
    aes(
      x = target_doublet_rate_pct,
      y = .data[[ycol]],
      color = method,
      group = method
    )
  ) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.4) +
    scale_color_manual(
      values = method_colors,
      labels = method_labels,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(5, 10, 20, 30, 40),
      labels = paste0(c(5, 10, 20, 30, 40), "%")
    ) +
    labs(
      title = title_text,
      x = "Doublet rate",
      y = ylab,
      color = "Method"
    ) +
    base_theme
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = number_format(accuracy = 0.01)
    )
  }
  
  if (show_legend) {
    p <- p + legend_inside_right
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

############################
## 6. Five line plots
############################

## (1) AUPRC
p_auprc <- make_line_plot(
  df = metrics_all,
  ycol = "AUPRC",
  ylab = "AUPRC",
  title_text = "AUPRC",
  show_legend = FALSE,
  y_limits = c(0.5, 1),
  y_breaks = seq(0.5, 1, 0.1)
)

## (2) AUROC
p_auroc <- make_line_plot(
  df = metrics_all,
  ycol = "AUROC",
  ylab = "AUROC",
  title_text = "AUROC",
  show_legend = TRUE,
  y_limits = c(0.85, 1),
  y_breaks = seq(0.85, 1, 0.05)
)

## (3) run_time
runtime_df <- metrics_all %>%
  filter(method %in% c("scDblFinder_random", "scds_cxds", "DoubletFinder"))
runtime_colors <- c(
  "scDblFinder_random" = "#E31A1B",
  "scds_cxds"          = "#A66527",  # scds 统一颜色
  "DoubletFinder"      = "#FD7F00"
)
p_runtime <- ggplot(
  runtime_df,
  aes(
    x = target_doublet_rate_pct,
    y = run_time,
    color = method,
    group = method
  )
) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.4) +
  scale_color_manual(
    values = runtime_colors,
    labels = c("scDblFinder_random" = "scDblFinder",
             "scds_cxds" = "scds",
             "DoubletFinder" = "DoubletFinder"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(5, 10, 20, 30, 40),
    labels = paste0(c(5, 10, 20, 30, 40), "%")
  ) +
  scale_y_continuous(
    limits = c(0,30),
    breaks = seq(0,30,5)
  )+
  labs(
    title = "Run time",
    x = "Target doublet rate",
    y = "Run time (sec)",
    color = "Method"
  ) +
  base_theme +
  legend_inside_right

## (4) homotypic recall
p_homo <- make_line_plot(
  df = metrics_all,
  ycol = "homotypic_recall",
  ylab = "Recall",
  title_text = "Homotypic recall",
  show_legend = FALSE,
  y_limits = c(0.2, 1),
  y_breaks = seq(0.2, 1, 0.2)
)

## (5) heterotypic recall
p_hetero <- make_line_plot(
  df = metrics_all,
  ycol = "heterotypic_recall",
  ylab = "Recall",
  title_text = "Heterotypic recall",
  show_legend = TRUE,
  y_limits = c(0.5, 1),
  y_breaks = seq(0.5, 1, 0.1)
)

############################
## 7. Heatmap: outperforming counts
##
## Count how many of the 5 datasets each method
## is the best-performing method on for:
## - AUPRC
## - AUROC
## - run_time (smaller is better)
## - homotypic recall
## - heterotypic recall
############################

metric_def <- tibble::tribble(
  ~metric_display,        ~metric_col,           ~direction,
  "AUPRC",                "AUPRC",               "max",
  "AUROC",                "AUROC",               "max",
  "Running time",             "run_time",            "min",
  "Homotypic recall",     "homotypic_recall",    "max",
  "Heterotypic recall",   "heterotypic_recall",  "max"
)

count_best_methods <- function(df, metric_col, direction = c("max", "min")) {
  direction <- match.arg(direction)
  
  tmp <- df %>%
    select(scenario, method, value = all_of(metric_col)) %>%
    filter(!is.na(value))
  
  tmp_best <- tmp %>%
    group_by(scenario) %>%
    mutate(
      best_value = if (direction == "max") max(value, na.rm = TRUE) else min(value, na.rm = TRUE),
      is_best = abs(value - best_value) < 1e-12
    ) %>%
    ungroup() %>%
    filter(is_best) %>%
    count(method, name = "n_best") %>%
    complete(
      method = factor(method_order, levels = method_order),
      fill = list(n_best = 0)
    ) %>%
    mutate(method = factor(method, levels = method_order))
  
  return(tmp_best)
}

heatmap_long <- lapply(seq_len(nrow(metric_def)), function(i) {
  row_i <- metric_def[i, ]
  
  count_best_methods(
    df = metrics_all,
    metric_col = row_i$metric_col,
    direction = row_i$direction
  ) %>%
    mutate(metric = row_i$metric_display)
}) %>%
  bind_rows() %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "AUPRC",
        "AUROC",
        "Running time",
        "Homotypic recall",
        "Heterotypic recall"
      )
    ),
    method = factor(method, levels = method_order),
    method_label = factor(method_labels[as.character(method)],
                          levels = method_labels[method_order])
  )

write.csv(
  heatmap_long,
  file.path(plot_dir, "outperforming_count_heatmap_data.csv"),
  row.names = FALSE
)

p_heatmap <- ggplot(
  heatmap_long,
  aes(
    x = method_label,
    y = metric,
    fill = n_best
  )
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = n_best), size = 4.2, fontface = "bold") +
  scale_fill_gradient(
    low = "#63be7b",
    high = "#FA8170",
    limits = c(0, 5),
    breaks = 0:5
  ) +
  labs(
    title = "Best-performing method counts",
    x = "Method",
    y = NULL,
    fill = "# best\n(out of 5)"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

############################
## 8. Save individual plots
############################
ggsave(file.path(plot_dir, "Fig1_AUPRC.pdf"), p_auprc, width = 8, height = 6, unit = "cm")

ggsave(file.path(plot_dir, "Fig1_AUROC.pdf"), p_auroc, width = 10.5, height = 6, unit = "cm")

ggsave(file.path(plot_dir, "Fig1_runtime.pdf"), p_runtime, width = 10.5, height = 6, unit = "cm")

ggsave(file.path(plot_dir, "Fig1_homotypic_recall.pdf"), p_homo, width = 8, height = 6, unit = "cm")

ggsave(file.path(plot_dir, "Fig1_heterotypic_recall.pdf"), p_hetero, width = 10.5, height = 6, unit="cm")

ggsave(file.path(plot_dir, "Fig1_outperforming_heatmap.pdf"), p_heatmap, width = 10, height = 7, unit = "cm")

############################
## 9. Combined 2 x 3 layout
##
## Row 1: AUPRC | AUROC | run_time
## Row 2: homotypic | heterotypic | heatmap
############################

row1 <- p_auprc + p_auroc + p_runtime + plot_layout(ncol = 3, widths = c(1, 1, 1))
row2 <- p_homo + p_hetero + p_heatmap + plot_layout(ncol = 3, widths = c(1, 1, 1))

p_combined <- row1 / row2 +
  plot_layout(heights = c(1, 1))

print(p_combined)

ggsave(
  filename = file.path(plot_dir, "Figure1_combined_2x3.png"),
  plot = p_combined,
  width = 14,
  height = 8.5,
  dpi = 300
)

ggsave(
  filename = file.path(plot_dir, "Figure1_combined_2x3.pdf"),
  plot = p_combined,
  width = 14,
  height = 8.5
)

############################
## 10. Also save combined metrics table
############################
write.csv(
  metrics_all,
  file.path(plot_dir, "combined_metrics_all_for_plotting.csv"),
  row.names = FALSE
)

cat("\nAll plots were saved to:\n")
cat(plot_dir, "\n")