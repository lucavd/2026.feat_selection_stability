# ==============================================================================
# Script 06 — Compute evaluation metrics
# ==============================================================================
# Purpose: Compute stability, accuracy, prediction, and parsimony metrics
#          from the feature selection results.
#
# Input:   results/feature_selection/fs_*.rds from Script 05
# Output:  results/metrics/metrics_all.rds (master table)
#          results/metrics/metrics_by_scenario.rds
#          results/metrics/metrics_summary.rds
# ==============================================================================

cat("=== Script 06: Metrics ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(data.table)
  library(pROC)
  library(future)
  library(furrr)
})

source(here("R", "utils", "helpers.R"))
source(here("R", "utils", "stability_metrics.R"))
config <- load_config()
setup_logging("Script 06: Metrics Computation", config)

fs_dir      <- get_path(config, "results_fs")
sim_dir     <- get_path(config, "simulated_data")
metrics_dir <- get_path(config, "results_metrics")

# --- List FS result files -----------------------------------------------------
fs_files <- list.files(fs_dir, pattern = "^fs_.*\\.rds$", full.names = TRUE)
cli::cli_alert_info("FS result files: {length(fs_files)}")

# ==============================================================================
# Metric computation functions
# ==============================================================================

#' Compute accuracy metrics (TPR, FDR, F1, MCC)
#' given selected features and ground truth
#'
#' @param selected Logical vector (p) of selected features
#' @param true_features Logical vector (p) of truly differential features
#' @return Named list of metrics
compute_accuracy_metrics <- function(selected, true_features) {
  tp <- sum(selected & true_features)
  fp <- sum(selected & !true_features)
  fn <- sum(!selected & true_features)
  tn <- sum(!selected & !true_features)

  tpr <- if (tp + fn > 0) tp / (tp + fn) else 0
  fdr <- if (tp + fp > 0) fp / (tp + fp) else 0
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  f1 <- if (precision + tpr > 0) 2 * precision * tpr / (precision + tpr) else 0

  # Matthews Correlation Coefficient
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (denom > 0) (tp * tn - fp * fn) / denom else 0

  list(tpr = tpr, fdr = fdr, precision = precision, f1 = f1, mcc = mcc,
       tp = tp, fp = fp, fn = fn, tn = tn)
}

#' Compute prediction AUC using selected features
#'
#' @param X Matrix (n × p)
#' @param y Binary response
#' @param selected Logical vector of selected features
#' @param test_fraction Fraction for test set
#' @param seed Random seed
#' @return AUC value
compute_prediction_metrics <- function(X, y, selected, test_fraction = 0.3, seed = 42) {
  result <- list(auc = NA_real_, balanced_accuracy = NA_real_)
  if (sum(selected) == 0) return(result)
  if (sum(selected) > nrow(X) - 5) {
    # More features than useful for logistic regression
    # Use only top features
    selected[selected] <- FALSE
    selected[seq_len(min(nrow(X) %/% 3, sum(selected)))] <- TRUE
    if (sum(selected) == 0) return(result)
  }

  set.seed(seed)
  n <- nrow(X)
  test_idx <- sample(n, round(n * test_fraction))
  train_idx <- setdiff(seq_len(n), test_idx)

  X_sel <- X[, selected, drop = FALSE]

  tryCatch({
    # Simple logistic regression
    df_train <- data.frame(y = y[train_idx], X_sel[train_idx, , drop = FALSE])
    df_test  <- data.frame(X_sel[test_idx, , drop = FALSE])

    fit <- glm(y ~ ., data = df_train, family = "binomial",
               control = list(maxit = 100))
    pred <- predict(fit, newdata = df_test, type = "response")
    roc_obj <- pROC::roc(y[test_idx], pred, quiet = TRUE)
    result$auc <- as.numeric(pROC::auc(roc_obj))

    # Balanced accuracy at optimal threshold
    pred_class <- as.integer(pred >= 0.5)
    y_test <- y[test_idx]
    sens <- sum(pred_class == 1 & y_test == 1) / max(sum(y_test == 1), 1)
    spec <- sum(pred_class == 0 & y_test == 0) / max(sum(y_test == 0), 1)
    result$balanced_accuracy <- (sens + spec) / 2

    result
  }, error = function(e) {
    result
  })
}

# ==============================================================================
# Process all FS results
# ==============================================================================

cli::cli_h1("Computing metrics")

setup_parallel(config)

# Group FS files by task_id (dataset) to aggregate across methods
# Parse file names: fs_<scenario>_<level>_rep<r>_<method>.rds
parse_fs_filename <- function(fname) {
  base <- gsub("^fs_|\\.rds$", "", basename(fname))
  # Last segment after _ is method name
  # But method names can contain underscores (e.g., "wilcoxon_fdr")
  # Load file to get metadata instead
  NULL
}

all_metrics <- furrr::future_map(fs_files, function(fpath) {
  tryCatch({
    fs_result <- readRDS(fpath)

    # Extract identifiers
    task_id     <- fs_result$task_id
    scenario    <- fs_result$scenario
    level_label <- fs_result$level_label
    rep_num     <- fs_result$rep
    method      <- fs_result$method
    category    <- fs_result$category
    true_feat   <- fs_result$true_features

    if (is.null(fs_result$selection_matrix)) {
      return(NULL)
    }

    sel_mat <- fs_result$selection_matrix
    imp_mat <- fs_result$importance_matrix

    # Only use converged runs
    if (!is.null(fs_result$converged)) {
      converged_idx <- which(fs_result$converged)
      if (length(converged_idx) < 2) {
        return(data.frame(
          task_id = task_id, scenario = scenario, level = level_label,
          rep = rep_num, method = method, category = category,
          status = "insufficient_convergence",
          stringsAsFactors = FALSE
        ))
      }
      sel_mat <- sel_mat[converged_idx, , drop = FALSE]
      if (!is.null(imp_mat)) imp_mat <- imp_mat[converged_idx, , drop = FALSE]
    }

    # --- Stability metrics ---------------------------------------------------
    stab <- compute_all_stability(sel_mat, imp_mat)

    # --- Accuracy metrics (average over bootstrap runs) ----------------------
    # Use the "stably selected" features (selected in >50% of bootstraps)
    selection_freq <- colMeans(sel_mat)
    for (threshold in config$metrics$stably_selected_thresholds) {
      stably_selected <- selection_freq >= threshold
      acc <- compute_accuracy_metrics(stably_selected, true_feat)
      stab[[paste0("tpr_t", threshold * 100)]] <- acc$tpr
      stab[[paste0("fdr_t", threshold * 100)]] <- acc$fdr
      stab[[paste0("f1_t", threshold * 100)]]  <- acc$f1
      stab[[paste0("mcc_t", threshold * 100)]] <- acc$mcc
      stab[[paste0("n_selected_t", threshold * 100)]] <- sum(stably_selected)
    }

    # Primary accuracy at 50% threshold
    stably_50 <- selection_freq >= 0.50
    acc_primary <- compute_accuracy_metrics(stably_50, true_feat)

    # --- Prediction metrics (using stably selected features) ------------------
    # Load simulated data for prediction
    sim_file <- file.path(sim_dir, paste0("sim_", task_id, ".rds"))
    pred_metrics <- list(auc = NA_real_, balanced_accuracy = NA_real_)
    if (file.exists(sim_file) && sum(stably_50) > 0) {
      sim_data <- readRDS(sim_file)
      pred_metrics <- compute_prediction_metrics(
        sim_data$X, sim_data$y, stably_50,
        test_fraction = config$metrics$prediction$test_fraction,
        seed = derive_seed(config$project$seed, paste0("auc_", task_id, "_", method))
      )
    }

    # --- Build result row ----------------------------------------------------
    data.frame(
      task_id = task_id,
      scenario = scenario,
      level = level_label,
      rep = rep_num,
      method = method,
      category = category,
      status = "success",
      # Stability
      nogueira = stab$nogueira$estimate,
      nogueira_ci_lower = stab$nogueira$ci_lower,
      nogueira_ci_upper = stab$nogueira$ci_upper,
      jaccard = stab$jaccard,
      kuncheva = stab$kuncheva,
      dice = stab$dice,
      spearman_rank = ifelse(is.null(stab$spearman_rank), NA_real_, stab$spearman_rank),
      # Accuracy (at 50% threshold)
      tpr = acc_primary$tpr,
      fdr = acc_primary$fdr,
      f1 = acc_primary$f1,
      mcc = acc_primary$mcc,
      # Prediction
      auc = pred_metrics$auc,
      balanced_accuracy = pred_metrics$balanced_accuracy,
      # Parsimony
      mean_n_selected = stab$mean_n_selected,
      sd_n_selected = stab$sd_n_selected,
      n_stably_selected_50 = sum(stably_50),
      # Bootstrap info
      n_bootstrap = stab$n_bootstrap,
      n_converged = fs_result$n_converged,
      total_time = fs_result$total_time,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cli::cli_alert_danger("Error processing {basename(fpath)}: {e$message}")
    NULL
  })
}, .options = furrr::furrr_options(seed = TRUE),
   .progress = TRUE)

teardown_parallel()

# ==============================================================================
# Combine and save
# ==============================================================================

cli::cli_h1("Aggregating results")

# Remove NULLs and bind
all_metrics <- all_metrics[!vapply(all_metrics, is.null, logical(1))]
metrics_dt <- data.table::rbindlist(all_metrics, fill = TRUE)

cli::cli_alert_info("Total metric rows: {nrow(metrics_dt)}")
cli::cli_alert_info("Successful: {sum(metrics_dt$status == 'success', na.rm = TRUE)}")

# Save master table
save_result(metrics_dt,
            file.path(metrics_dir, "metrics_all.rds"),
            metadata = list(script = "06_metrics.R",
                            n_rows = nrow(metrics_dt)))

# --- Aggregated summaries per scenario × method ------------------------------
if (nrow(metrics_dt[metrics_dt$status == "success", ]) > 0) {
  metrics_success <- metrics_dt[status == "success"]

  summary_dt <- metrics_success[, .(
    # Stability
    nogueira_mean = mean(nogueira, na.rm = TRUE),
    nogueira_sd   = sd(nogueira, na.rm = TRUE),
    jaccard_mean  = mean(jaccard, na.rm = TRUE),
    # Accuracy
    tpr_mean = mean(tpr, na.rm = TRUE),
    tpr_sd   = sd(tpr, na.rm = TRUE),
    fdr_mean = mean(fdr, na.rm = TRUE),
    fdr_sd   = sd(fdr, na.rm = TRUE),
    f1_mean  = mean(f1, na.rm = TRUE),
    mcc_mean = mean(mcc, na.rm = TRUE),
    # Prediction
    auc_mean = mean(auc, na.rm = TRUE),
    auc_sd   = sd(auc, na.rm = TRUE),
    # Parsimony
    n_selected_mean = mean(mean_n_selected, na.rm = TRUE),
    n_stably_50_mean = mean(n_stably_selected_50, na.rm = TRUE),
    # Meta
    n_reps = .N,
    total_time = sum(total_time, na.rm = TRUE)
  ), by = .(scenario, level, method, category)]

  save_result(summary_dt,
              file.path(metrics_dir, "metrics_summary.rds"),
              metadata = list(script = "06_metrics.R"))

  # Print summary table
  cli::cli_h2("Summary by method (averaged across all scenarios)")
  method_summary <- metrics_success[, .(
    nogueira = round(mean(nogueira, na.rm = TRUE), 3),
    jaccard = round(mean(jaccard, na.rm = TRUE), 3),
    tpr = round(mean(tpr, na.rm = TRUE), 3),
    fdr = round(mean(fdr, na.rm = TRUE), 3),
    auc = round(mean(auc, na.rm = TRUE), 3),
    n_selected = round(mean(mean_n_selected, na.rm = TRUE), 1)
  ), by = .(method, category)]

  print(method_summary[order(category, method)])
}

finalize_logging("Script 06", session_file = "results/session_info_06.txt")

cat("\n=== Script 06: Complete ===\n")
