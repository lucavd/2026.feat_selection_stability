# ==============================================================================
# Script 07 — Cross-database validation on real data
# ==============================================================================
# Purpose: Apply all 12 FS methods to the real processed datasets with
#          bootstrap resampling, and cross-validate findings between
#          datasets of the same disease/platform.
#
# This script is INDEPENDENT from the simulation pipeline (04-06).
#
# Input:   data/processed/<accession>_processed.rds from Script 02
# Output:  results/cross_validation/cv_<accession>_<method>.rds
#          results/cross_validation/cv_summary.rds
# ==============================================================================

cat("=== Script 07: Cross-Database Validation ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(data.table)
  library(future)
  library(furrr)
})

source(here("R", "utils", "helpers.R"))
source(here("R", "utils", "fs_methods.R"))
source(here("R", "utils", "stability_metrics.R"))
config <- load_config()
setup_logging("Script 07: Cross-Database Validation", config)

proc_dir <- get_path(config, "processed_data")
cv_dir   <- get_path(config, "results_cv")

# --- Load processed datasets --------------------------------------------------
summary <- load_result(file.path(proc_dir, "datasets_summary.rds"))$data
successful_ids <- names(summary)[vapply(summary, function(s) s$status == "success",
                                        logical(1))]

cli::cli_alert_info("{length(successful_ids)} datasets available for cross-validation")

methods     <- config$methods
n_bootstrap <- config$simulation$resampling$n_bootstrap
subsample_frac <- config$simulation$resampling$subsample_fraction

# ==============================================================================
# Part 1: Bootstrap FS on each real dataset
# ==============================================================================

cli::cli_h1("Part 1: Bootstrap FS on real datasets")

setup_parallel(config)

for (acc_id in successful_ids) {
  cli::cli_h2("Dataset: {acc_id}")

  dat <- load_result(file.path(proc_dir, paste0(acc_id, "_processed.rds")))$data
  X <- dat$X
  y <- dat$y

  for (method in methods) {
    cv_id <- paste0(acc_id, "_", method$name)
    out_file <- file.path(cv_dir, paste0("cv_", cv_id, ".rds"))

    if (file.exists(out_file)) {
      cli::cli_alert_info("  {method$name}: already computed, skipping")
      next
    }

    cli::cli_alert_info("  {method$name}: running {n_bootstrap} bootstraps...")

    base_seed <- derive_seed(config$project$seed, paste0("cv_", cv_id))

    # Run bootstrap FS
    p <- ncol(X)
    selection_matrix  <- matrix(0L, nrow = n_bootstrap, ncol = p)
    importance_matrix <- matrix(NA_real_, nrow = n_bootstrap, ncol = p)
    colnames(selection_matrix)  <- colnames(X)
    colnames(importance_matrix) <- colnames(X)
    timings   <- numeric(n_bootstrap)
    converged <- logical(n_bootstrap)

    for (b in seq_len(n_bootstrap)) {
      set.seed(base_seed + b)

      # Stratified subsample
      idx_0 <- which(y == 0)
      idx_1 <- which(y == 1)
      n_sub_0 <- max(round(length(idx_0) * subsample_frac), 2)
      n_sub_1 <- max(round(length(idx_1) * subsample_frac), 2)
      boot_idx <- c(
        sample(idx_0, min(n_sub_0, length(idx_0)), replace = FALSE),
        sample(idx_1, min(n_sub_1, length(idx_1)), replace = FALSE)
      )

      result <- run_fs_method(method$name, X[boot_idx, ], y[boot_idx], method$params)

      sel <- result$selected
      sel[is.na(sel)] <- FALSE
      selection_matrix[b, ]  <- as.integer(sel)
      importance_matrix[b, ] <- result$importance
      timings[b]   <- result$time
      converged[b] <- result$converged
    }

    n_converged <- sum(converged)
    if (n_converged >= 2) {
      stab <- compute_all_stability(
        selection_matrix[converged, , drop = FALSE],
        importance_matrix[converged, , drop = FALSE]
      )
      selection_freq <- colMeans(selection_matrix[converged, , drop = FALSE])
    } else if (n_converged == 1) {
      selection_freq <- as.numeric(selection_matrix[converged, , drop = FALSE])
      names(selection_freq) <- colnames(X)
      stab <- list(
        nogueira = list(estimate = NA_real_, variance = NA_real_,
                        ci_lower = NA_real_, ci_upper = NA_real_),
        jaccard = NA_real_,
        kuncheva = NA_real_,
        dice = NA_real_,
        spearman_rank = NA_real_,
        n_bootstrap = 1L,
        p = p,
        mean_n_selected = mean(rowSums(selection_matrix[converged, , drop = FALSE])),
        sd_n_selected = NA_real_
      )
    } else {
      selection_freq <- rep(NA_real_, p)
      names(selection_freq) <- colnames(X)
      stab <- list(
        nogueira = list(estimate = NA_real_, variance = NA_real_,
                        ci_lower = NA_real_, ci_upper = NA_real_),
        jaccard = NA_real_,
        kuncheva = NA_real_,
        dice = NA_real_,
        spearman_rank = NA_real_,
        n_bootstrap = 0L,
        p = p,
        mean_n_selected = NA_real_,
        sd_n_selected = NA_real_
      )
    }

    cv_result <- list(
      accession = acc_id,
      method = method$name,
      category = method$category,
      selection_matrix = selection_matrix,
      importance_matrix = importance_matrix,
      selection_freq = selection_freq,
      stability = stab,
      n_converged = n_converged,
      total_time = sum(timings),
      feature_names = colnames(X)
    )

    saveRDS(cv_result, out_file)
    cli::cli_alert_success("  {method$name}: nogueira={round(stab$nogueira$estimate, 3)}, n_selected_50={sum(selection_freq >= config$metrics$stably_selected_thresholds[1], na.rm = TRUE)}")
  }
}

# ==============================================================================
# Part 2: Cross-database concordance
# ==============================================================================

cli::cli_h1("Part 2: Cross-database concordance")

# Find pairs of datasets that can be compared
# (same disease/condition or same platform)
dataset_info <- lapply(successful_ids, function(id) {
  dat <- load_result(file.path(proc_dir, paste0(id, "_processed.rds")))$data
  list(id = id, platform = dat$platform, description = dat$description,
       features = colnames(dat$X))
})
names(dataset_info) <- successful_ids

# Build pairs by platform
platform_groups <- split(successful_ids,
                         vapply(dataset_info, function(d) d$platform, character(1)))

concordance_results <- list()

for (plat in names(platform_groups)) {
  ids <- platform_groups[[plat]]
  if (length(ids) < 2) next

  cli::cli_h2("Platform: {plat} ({length(ids)} datasets)")

  pairs <- combn(ids, 2)

  for (pair_idx in seq_len(ncol(pairs))) {
    id1 <- pairs[1, pair_idx]
    id2 <- pairs[2, pair_idx]
    cli::cli_alert_info("Comparing {id1} vs {id2}")

    for (method in methods) {
      cv_file1 <- file.path(cv_dir, paste0("cv_", id1, "_", method$name, ".rds"))
      cv_file2 <- file.path(cv_dir, paste0("cv_", id2, "_", method$name, ".rds"))

      if (!file.exists(cv_file1) || !file.exists(cv_file2)) next

      cv1 <- readRDS(cv_file1)
      cv2 <- readRDS(cv_file2)

      # Find common features (by name)
      common_features <- intersect(cv1$feature_names, cv2$feature_names)
      if (length(common_features) < 10) {
        cli::cli_alert_warning("  {method$name}: only {length(common_features)} common features, skipping")
        next
      }

      # Compare selection frequencies on common features
      freq1 <- cv1$selection_freq[common_features]
      freq2 <- cv2$selection_freq[common_features]

      # Concordance metrics
      # 1. Spearman correlation of selection frequencies
      spearman_cor <- cor(freq1, freq2, method = "spearman", use = "complete.obs")

      # 2. Jaccard of stably selected features (>50%)
      sel1 <- names(freq1)[which(!is.na(freq1) & freq1 >= config$metrics$stably_selected_thresholds[1])]
      sel2 <- names(freq2)[which(!is.na(freq2) & freq2 >= config$metrics$stably_selected_thresholds[1])]
      jaccard <- length(intersect(sel1, sel2)) /
        max(length(union(sel1, sel2)), 1)

      # 3. Overlap coefficient
      overlap <- length(intersect(sel1, sel2)) /
        max(min(length(sel1), length(sel2)), 1)

      pair_id <- paste0(id1, "_vs_", id2)
      concordance_results[[paste0(pair_id, "_", method$name)]] <- data.frame(
        dataset1 = id1,
        dataset2 = id2,
        platform = plat,
        method = method$name,
        category = method$category,
        n_common_features = length(common_features),
        spearman_cor = spearman_cor,
        jaccard = jaccard,
        overlap = overlap,
        n_selected_1 = length(sel1),
        n_selected_2 = length(sel2),
        n_overlap = length(intersect(sel1, sel2)),
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine concordance results
if (length(concordance_results) > 0) {
  concordance_dt <- data.table::rbindlist(concordance_results, fill = TRUE)
  cli::cli_alert_info("Concordance pairs computed: {nrow(concordance_dt)}")

  save_result(concordance_dt,
              file.path(cv_dir, "concordance_results.rds"),
              metadata = list(script = "07_cross_validation.R"))
}

# ==============================================================================
# Summary
# ==============================================================================

cli::cli_h1("Cross-Validation Summary")

# Per-method stability on real data
cv_files <- list.files(cv_dir, pattern = "^cv_.*\\.rds$", full.names = TRUE)
cv_files <- cv_files[!grepl("summary|concordance", cv_files)]

cv_summary_list <- lapply(cv_files, function(f) {
  tryCatch({
    res <- readRDS(f)
    data.frame(
      accession = res$accession,
      method = res$method,
      category = res$category,
      nogueira = res$stability$nogueira$estimate,
      jaccard = res$stability$jaccard,
      n_selected_50 = sum(res$selection_freq >= config$metrics$stably_selected_thresholds[1], na.rm = TRUE),
      n_converged = res$n_converged,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
})

cv_summary_dt <- data.table::rbindlist(
  cv_summary_list[!vapply(cv_summary_list, is.null, logical(1))],
  fill = TRUE
)

if (nrow(cv_summary_dt) > 0) {
  save_result(cv_summary_dt,
              file.path(cv_dir, "cv_summary.rds"),
              metadata = list(script = "07_cross_validation.R"))

  cli::cli_h2("Stability on real data by method")
  method_avg <- cv_summary_dt[, .(
    nogueira_mean = round(mean(nogueira, na.rm = TRUE), 3),
    jaccard_mean  = round(mean(jaccard, na.rm = TRUE), 3),
    n_selected    = round(mean(n_selected_50, na.rm = TRUE), 1)
  ), by = .(method, category)]

  print(method_avg[order(category, method)])
}

teardown_parallel()
finalize_logging("Script 07", session_file = "results/session_info_07.txt")

cat("\n=== Script 07: Complete ===\n")
