# ==============================================================================
# Script 03 — Extract empirical parameters from real datasets
# ==============================================================================
# Purpose: Estimate distributional and correlation parameters from processed
#          real-world metabolomics data to inform realistic simulations.
# Input:   data/processed/<accession>_processed.rds from Script 02
# Output:  data/empirical_params/empirical_params.rds
# ==============================================================================

cat("=== Script 03: Extract Empirical Parameters ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(data.table)
  library(mvtnorm)
  library(corpcor)
  library(MASS)
})

source(here("R", "utils", "helpers.R"))
config <- load_config()
setup_logging("Script 03: Extract Empirical Parameters", config)

proc_dir   <- get_path(config, "processed_data")
params_dir <- get_path(config, "empirical_params")

# --- Load processed datasets --------------------------------------------------
summary <- load_result(file.path(proc_dir, "datasets_summary.rds"))$data
successful_ids <- names(summary)[vapply(summary, function(s) s$status == "success",
                                        logical(1))]

cli::cli_alert_info("Extracting parameters from {length(successful_ids)} datasets")

if (length(successful_ids) == 0) {
  cli::cli_abort("No processed datasets available. Run Script 02 first.")
}

# ==============================================================================
# Parameter extraction per dataset
# ==============================================================================

all_params <- list()

for (acc_id in successful_ids) {
  cli::cli_h2("Extracting: {acc_id}")

  dat <- load_result(file.path(proc_dir, paste0(acc_id, "_processed.rds")))$data
  X <- dat$X
  y <- dat$y
  n <- nrow(X)
  p <- ncol(X)

  params <- list(
    id = acc_id,
    source = dat$source,
    platform = dat$platform,
    n = n,
    p = p,
    n_case = sum(y == 1),
    n_control = sum(y == 0)
  )

  # --- 1. Marginal distributions -------------------------------------------
  # Mean and SD per feature, stratified by group
  params$mean_control <- colMeans(X[y == 0, , drop = FALSE], na.rm = TRUE)
  params$mean_case    <- colMeans(X[y == 1, , drop = FALSE], na.rm = TRUE)
  params$sd_control   <- apply(X[y == 0, , drop = FALSE], 2, sd, na.rm = TRUE)
  params$sd_case      <- apply(X[y == 1, , drop = FALSE], 2, sd, na.rm = TRUE)

  # Overall
  params$mean_overall <- colMeans(X, na.rm = TRUE)
  params$sd_overall   <- apply(X, 2, sd, na.rm = TRUE)

  # Observed fold changes (data is on log scale, so FC = difference, not ratio)
  params$observed_log2fc <- (params$mean_case - params$mean_control) / log(2)

  cli::cli_alert_info("  Median |log2FC|: {round(median(abs(params$observed_log2fc), na.rm=TRUE), 3)}")

  # --- 2. Correlation structure ---------------------------------------------
  # Estimate correlation matrix from control group (unperturbed)
  X_ctrl <- X[y == 0, , drop = FALSE]

  # Remove features with zero variance or all-NA in controls (unusable for cor)
  ctrl_sds <- apply(X_ctrl, 2, sd, na.rm = TRUE)
  bad_cols <- is.na(ctrl_sds) | ctrl_sds < 1e-10
  params$cor_features_removed <- sum(bad_cols)
  params$cor_features_kept <- which(!bad_cols)  # Indices into original feature space
  if (any(bad_cols)) {
    cli::cli_alert_warning("  Removing {sum(bad_cols)} zero-variance/all-NA features from correlation estimation")
    X_ctrl <- X_ctrl[, !bad_cols, drop = FALSE]
  }
  n_ctrl <- nrow(X_ctrl)
  p_ctrl <- ncol(X_ctrl)

  # Use shrinkage estimator for stability (Ledoit-Wolf via corpcor)
  if (p_ctrl > n_ctrl) {
    # High-dimensional: shrinkage is critical
    cli::cli_alert_info("  p > n: using shrinkage correlation estimator")
    cor_shrink <- corpcor::cor.shrink(X_ctrl, verbose = FALSE)
    params$cor_matrix <- as.matrix(cor_shrink)
    params$cor_method <- "ledoit_wolf_shrinkage"
  } else {
    # Standard correlation with shrinkage as safety net
    cor_raw <- cor(X_ctrl, use = "pairwise.complete.obs")
    # Check for NA/NaN in correlation matrix (e.g. from remaining missing data)
    if (any(is.na(cor_raw))) {
      cli::cli_alert_info("  NA values in correlation matrix; falling back to shrinkage")
      cor_shrink <- corpcor::cor.shrink(X_ctrl, verbose = FALSE)
      params$cor_matrix <- as.matrix(cor_shrink)
      params$cor_method <- "ledoit_wolf_shrinkage"
    } else {
      # Check if positive definite
      eigenvals <- eigen(cor_raw, symmetric = TRUE, only.values = TRUE)$values
      if (any(eigenvals <= 0)) {
        cli::cli_alert_info("  Correlation matrix not PD; applying shrinkage")
        cor_shrink <- corpcor::cor.shrink(X_ctrl, verbose = FALSE)
        params$cor_matrix <- as.matrix(cor_shrink)
        params$cor_method <- "ledoit_wolf_shrinkage"
      } else {
        params$cor_matrix <- cor_raw
        params$cor_method <- "pearson"
      }
    }
  }
  params$n_features_cor <- p_ctrl  # Track how many features used for correlation

  # Summary statistics for correlation structure
  cor_vals <- params$cor_matrix[lower.tri(params$cor_matrix)]
  params$cor_summary <- list(
    mean_abs_cor   = mean(abs(cor_vals)),
    median_abs_cor = median(abs(cor_vals)),
    q25_abs_cor    = quantile(abs(cor_vals), 0.25),
    q75_abs_cor    = quantile(abs(cor_vals), 0.75),
    max_abs_cor    = max(abs(cor_vals)),
    prop_above_03  = mean(abs(cor_vals) > 0.3),
    prop_above_05  = mean(abs(cor_vals) > 0.5)
  )

  cli::cli_alert_info("  Mean |cor|: {round(params$cor_summary$mean_abs_cor, 3)}")
  cli::cli_alert_info("  % |cor| > 0.3: {round(params$cor_summary$prop_above_03 * 100, 1)}%")

  # --- 3. Eigenvalue spectrum -----------------------------------------------
  # Captures the effective dimensionality / correlation decay
  eigenvals <- eigen(params$cor_matrix, symmetric = TRUE, only.values = TRUE)$values
  eigenvals <- pmax(eigenvals, 0)  # Ensure non-negative
  params$eigenvalues <- eigenvals
  params$effective_rank <- sum(eigenvals > 1)  # Kaiser criterion
  params$explained_var_top10 <- sum(eigenvals[1:min(10, length(eigenvals))]) / sum(eigenvals)

  cli::cli_alert_info("  Effective rank: {params$effective_rank}")
  cli::cli_alert_info("  Top 10 eigenvalues explain: {round(params$explained_var_top10 * 100, 1)}%")

  # --- 4. Skewness and kurtosis (per feature) --------------------------------
  params$skewness <- apply(X, 2, function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 3) return(NA_real_)
    m <- mean(x)
    s <- sd(x)
    if (is.na(s) || s < 1e-10) return(0)
    (n / ((n-1)*(n-2))) * sum(((x - m) / s)^3)
  })

  params$kurtosis <- apply(X, 2, function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 4) return(NA_real_)
    m <- mean(x)
    s <- sd(x)
    if (is.na(s) || s < 1e-10) return(0)
    ((n*(n+1)) / ((n-1)*(n-2)*(n-3))) * sum(((x - m) / s)^4) -
      (3*(n-1)^2) / ((n-2)*(n-3))
  })

  cli::cli_alert_info("  Median skewness: {round(median(params$skewness, na.rm = TRUE), 3)}")

  # --- 5. Missing data patterns (from raw, before imputation) ----------------
  X_raw <- dat$X_raw
  params$missing_rate_overall <- mean(is.na(X_raw))
  params$missing_rate_per_feature <- colMeans(is.na(X_raw))

  cli::cli_alert_info("  Missing rate: {round(params$missing_rate_overall * 100, 1)}%")

  all_params[[acc_id]] <- params
  cli::cli_alert_success("  {acc_id}: parameters extracted")
}

# ==============================================================================
# Aggregate empirical parameters (pooled across datasets)
# ==============================================================================

cli::cli_h1("Aggregating parameters")

# Pool correlation matrices by platform
platforms <- unique(vapply(all_params, function(p) p$platform, character(1)))
pooled <- list()

for (plat in platforms) {
  plat_ids <- names(all_params)[vapply(all_params, function(p) p$platform == plat, logical(1))]
  cli::cli_alert_info("Platform {plat}: {length(plat_ids)} datasets")

  # Use the largest dataset's correlation as reference
  sizes <- vapply(plat_ids, function(id) all_params[[id]]$p, integer(1))
  ref_id <- plat_ids[which.max(sizes)]
  cli::cli_alert_info("  Reference dataset: {ref_id} (p={max(sizes)})")

  pooled[[plat]] <- list(
    reference_id = ref_id,
    cor_matrix = all_params[[ref_id]]$cor_matrix,
    mean_overall = all_params[[ref_id]]$mean_overall,
    sd_overall = all_params[[ref_id]]$sd_overall,
    eigenvalues = all_params[[ref_id]]$eigenvalues,
    n_datasets = length(plat_ids),
    dataset_ids = plat_ids
  )
}

# ==============================================================================
# Save
# ==============================================================================

results <- list(
  per_dataset = all_params,
  pooled = pooled,
  platforms = platforms
)

save_result(results,
            file.path(params_dir, "empirical_params.rds"),
            metadata = list(script = "03_extract_empirical_params.R",
                            n_datasets = length(all_params),
                            platforms = platforms))

finalize_logging("Script 03", session_file = "results/session_info_03.txt")

cat("\n=== Script 03: Complete ===\n")
