# ==============================================================================
# Script 08 — Figures and tables for the paper
# ==============================================================================
# Purpose: Generate publication-quality figures from simulation and
#          cross-validation results.
#
# Input:   results/metrics/metrics_all.rds (Script 06)
#          results/metrics/metrics_summary.rds (Script 06)
#          results/cross_validation/cv_summary.rds (Script 07)
#          results/cross_validation/concordance_results.rds (Script 07)
# Output:  results/figures/fig_*.pdf
# ==============================================================================

cat("=== Script 08: Figures ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
  library(ggrepel)
})

source(here("R", "utils", "helpers.R"))
config <- load_config()
setup_logging("Script 08: Figures", config)

metrics_dir <- get_path(config, "results_metrics")
cv_dir      <- get_path(config, "results_cv")
fig_dir     <- get_path(config, "figures")

# --- Figure settings ----------------------------------------------------------
fig_cfg <- config$figures
theme_paper <- theme_bw(base_size = fig_cfg$font_size) +
  theme(
    text = element_text(family = fig_cfg$font_family),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )
theme_set(theme_paper)

# Method category colors
cat_colors <- unlist(fig_cfg$method_colors)

# Helper to save figure
save_fig <- function(p, name, width = fig_cfg$width, height = fig_cfg$height) {
  fpath <- file.path(fig_dir, paste0(name, ".", fig_cfg$format))
  ggsave(fpath, plot = p, width = width, height = height, dpi = fig_cfg$dpi)
  cli::cli_alert_success("Saved: {name}.{fig_cfg$format}")
}

# --- Load data ----------------------------------------------------------------
metrics_all <- tryCatch(
  load_result(file.path(metrics_dir, "metrics_all.rds"))$data,
  error = function(e) {
    cli::cli_alert_danger("metrics_all.rds not found: {e$message}")
    NULL
  }
)

metrics_summary <- tryCatch(
  load_result(file.path(metrics_dir, "metrics_summary.rds"))$data,
  error = function(e) NULL
)

cv_summary <- tryCatch(
  load_result(file.path(cv_dir, "cv_summary.rds"))$data,
  error = function(e) NULL
)

concordance <- tryCatch(
  load_result(file.path(cv_dir, "concordance_results.rds"))$data,
  error = function(e) NULL
)

if (is.null(metrics_all)) {
  cli::cli_abort("No metrics data available. Run Scripts 05-06 first.")
}

metrics_ok <- metrics_all[status == "success"]

# ==============================================================================
# Figure 1: Stability overview — Nogueira index by method across all scenarios
# ==============================================================================

cli::cli_h1("Figure 1: Stability Overview")

fig1 <- ggplot(metrics_ok, aes(x = reorder(method, nogueira, FUN = median),
                                y = nogueira, fill = category)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = cat_colors, name = "Category") +
  coord_flip() +
  labs(x = NULL, y = "Nogueira Stability Index",
       title = "Feature selection stability across simulation scenarios") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.5, y = 0.52, label = "chance level",
           hjust = 0, size = 3, color = "grey50")

save_fig(fig1, "fig01_stability_overview", width = 8, height = 6)

# ==============================================================================
# Figure 2: Stability-Accuracy trade-off
# ==============================================================================

cli::cli_h1("Figure 2: Stability vs Accuracy")

if (!is.null(metrics_summary)) {
  fig2_data <- metrics_summary[, .(
    nogueira = mean(nogueira_mean, na.rm = TRUE),
    tpr = mean(tpr_mean, na.rm = TRUE),
    fdr = mean(fdr_mean, na.rm = TRUE),
    n_selected = mean(n_selected_mean, na.rm = TRUE)
  ), by = .(method, category)]

  fig2 <- ggplot(fig2_data, aes(x = nogueira, y = tpr, color = category)) +
    geom_point(aes(size = n_selected), alpha = 0.8) +
    geom_text_repel(aes(label = method), size = 3, max.overlaps = 15) +
    scale_color_manual(values = cat_colors, name = "Category") +
    scale_size_continuous(name = "Avg. features\nselected", range = c(2, 8)) +
    labs(x = "Stability (Nogueira Index)", y = "Sensitivity (TPR)",
         title = "Stability–accuracy trade-off") +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3)

  save_fig(fig2, "fig02_stability_accuracy_tradeoff", width = 8, height = 6)
}

# ==============================================================================
# Figure 3: Scenario-specific stability (heatmap-style faceted plot)
# ==============================================================================

cli::cli_h1("Figure 3: Scenario × Method Heatmap")

if (!is.null(metrics_summary)) {
  fig3 <- ggplot(metrics_summary, aes(x = method, y = level, fill = nogueira_mean)) +
    geom_tile(color = "white") +
    scale_fill_viridis(name = "Nogueira\nIndex", option = "D", limits = c(0, 1)) +
    facet_wrap(~scenario, scales = "free_y", ncol = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    labs(x = NULL, y = NULL,
         title = "Stability across simulation scenarios")

  save_fig(fig3, "fig03_scenario_heatmap", width = 12, height = 10)
}

# ==============================================================================
# Figure 4: Effect of p/n ratio on stability
# ==============================================================================

cli::cli_h1("Figure 4: p/n Ratio Effect")

s1_data <- metrics_ok[scenario == "S1_pn_ratio"]
if (nrow(s1_data) > 0) {
  fig4 <- ggplot(s1_data, aes(x = level, y = nogueira, fill = category)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = cat_colors, name = "Category") +
    facet_wrap(~method, ncol = 4) +
    labs(x = "p/n ratio", y = "Nogueira Stability Index",
         title = "Impact of dimensionality on feature selection stability") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_fig(fig4, "fig04_pn_ratio", width = 12, height = 8)
}

# ==============================================================================
# Figure 5: Effect of correlation structure
# ==============================================================================

cli::cli_h1("Figure 5: Correlation Effect")

s3_data <- metrics_ok[scenario == "S3_correlation"]
if (nrow(s3_data) > 0) {
  fig5 <- ggplot(s3_data, aes(x = level, y = nogueira, color = method, group = method)) +
    stat_summary(fun = mean, geom = "line", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, alpha = 0.5) +
    scale_color_viridis_d(name = "Method") +
    labs(x = "Correlation Level", y = "Nogueira Stability Index",
         title = "Impact of feature correlation on stability")

  save_fig(fig5, "fig05_correlation_effect", width = 9, height = 6)
}

# ==============================================================================
# Figure 6: FDR control
# ==============================================================================

cli::cli_h1("Figure 6: FDR Control")

fig6 <- ggplot(metrics_ok, aes(x = reorder(method, fdr, FUN = median),
                                y = fdr, fill = category)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = cat_colors, name = "Category") +
  coord_flip() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 0.5, y = 0.07, label = "nominal FDR = 0.05",
           hjust = 0, size = 3, color = "red") +
  labs(x = NULL, y = "False Discovery Rate",
       title = "FDR control across methods")

save_fig(fig6, "fig06_fdr_control", width = 8, height = 6)

# ==============================================================================
# Figure 7: Cross-database validation (real data)
# ==============================================================================

cli::cli_h1("Figure 7: Real Data Stability")

if (!is.null(cv_summary) && nrow(cv_summary) > 0) {
  fig7 <- ggplot(cv_summary, aes(x = reorder(method, nogueira, FUN = median),
                                  y = nogueira, fill = category)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = cat_colors, name = "Category") +
    coord_flip() +
    labs(x = NULL, y = "Nogueira Stability Index",
         title = "Feature selection stability on real metabolomics data")

  save_fig(fig7, "fig07_real_data_stability", width = 8, height = 6)
}

# ==============================================================================
# Figure 8: Cross-database concordance
# ==============================================================================

cli::cli_h1("Figure 8: Cross-Database Concordance")

if (!is.null(concordance) && nrow(concordance) > 0) {
  fig8 <- ggplot(concordance, aes(x = reorder(method, spearman_cor, FUN = median),
                                   y = spearman_cor, fill = category)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = cat_colors, name = "Category") +
    coord_flip() +
    labs(x = NULL, y = "Spearman Correlation of Selection Frequencies",
         title = "Cross-database concordance of feature selection") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3)

  save_fig(fig8, "fig08_cross_database_concordance", width = 8, height = 6)
}

# ==============================================================================
# Supplementary: Multi-panel composite
# ==============================================================================

cli::cli_h1("Supplementary Figure: Composite")

if (exists("fig1") && exists("fig6")) {
  fig_supp <- fig1 + fig6 +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Feature Selection Stability and FDR Control",
      tag_levels = "A"
    )
  save_fig(fig_supp, "fig_supp_composite", width = 14, height = 6)
}

# ==============================================================================
# Summary
# ==============================================================================

cli::cli_h1("Figures Summary")

fig_files <- list.files(fig_dir, pattern = paste0("\\.", fig_cfg$format, "$"))
cli::cli_alert_success("Generated {length(fig_files)} figures in {fig_dir}")
for (f in fig_files) {
  cli::cli_alert_info("  {f}")
}

finalize_logging("Script 08", session_file = "results/session_info_08.txt")

cat("\n=== Script 08: Complete ===\n")
