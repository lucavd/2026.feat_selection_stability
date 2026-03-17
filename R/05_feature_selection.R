# ==============================================================================
# Script 05 — Feature selection pipeline (BOTTLENECK)
# ==============================================================================
# Purpose: Apply all 12 FS methods to each simulated dataset with bootstrap
#          resampling to measure selection stability.
#
# Design:  For each simulated dataset:
#            For each FS method:
#              For each bootstrap sample (B=100):
#                - Draw stratified subsample (63.2%)
#                - Run FS method
#                - Record selected features + importance scores
#              Save selection_matrix (B × p) per method
#
# Estimated load: 7 scenarios × ~4 levels × 30 reps × 12 methods × 100 boot
#                 = ~1,008,000 fits → checkpoint aggressively
#
# Input:   data/simulated/sim_*.rds from Script 04
# Output:  results/feature_selection/fs_<task_id>_<method>.rds
# ==============================================================================

cat("=== Script 05: Feature Selection Pipeline ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(future)
  library(furrr)
  library(progressr)
})

source(here("R", "utils", "helpers.R"))
source(here("R", "utils", "fs_methods.R"))
config <- load_config()
setup_logging("Script 05: Feature Selection (BOTTLENECK)", config)

sim_dir <- get_path(config, "simulated_data")
fs_dir  <- get_path(config, "results_fs")

# --- Configuration ------------------------------------------------------------
methods     <- config$methods
n_bootstrap <- config$simulation$resampling$n_bootstrap
subsample_frac <- config$simulation$resampling$subsample_fraction
stratified  <- config$simulation$resampling$stratified

# --- List simulated datasets --------------------------------------------------
sim_files <- list.files(sim_dir, pattern = "^sim_.*\\.rds$", full.names = TRUE)
cli::cli_alert_info("Simulated datasets: {length(sim_files)}")
cli::cli_alert_info("Methods: {length(methods)}")
cli::cli_alert_info("Bootstrap samples: {n_bootstrap}")
cli::cli_alert_info("Total fits expected: ~{length(sim_files) * length(methods) * n_bootstrap}")

# ==============================================================================
# Bootstrap feature selection for a single dataset × method
# ==============================================================================

#' Run bootstrap FS for one dataset and one method
#'
#' @param sim_data Simulated dataset list (X, y, true_features, ...)
#' @param method_config Method configuration from config$methods
#' @param n_boot Number of bootstrap samples
#' @param subsample_frac Fraction of samples per bootstrap
#' @param stratified Logical, stratified sampling
#' @param base_seed Base seed for reproducibility
#' @return List with selection_matrix, importance_matrix, timings, convergence
run_bootstrap_fs <- function(sim_data, method_config, n_boot, subsample_frac,
                             stratified, base_seed) {
  X <- sim_data$X
  y <- sim_data$y
  n <- nrow(X)
  p <- ncol(X)

  method_name <- method_config$name
  params <- method_config$params

  # Pre-allocate
  selection_matrix  <- matrix(0L, nrow = n_boot, ncol = p)
  importance_matrix <- matrix(NA_real_, nrow = n_boot, ncol = p)
  colnames(selection_matrix)  <- colnames(X)
  colnames(importance_matrix) <- colnames(X)
  timings     <- numeric(n_boot)
  converged   <- logical(n_boot)
  messages    <- character(n_boot)

  for (b in seq_len(n_boot)) {
    set.seed(base_seed + b)

    # --- Stratified subsample -----------------------------------------------
    n_sub <- round(n * subsample_frac)
    if (stratified) {
      idx_0 <- which(y == 0)
      idx_1 <- which(y == 1)
      n_sub_0 <- round(length(idx_0) * subsample_frac)
      n_sub_1 <- round(length(idx_1) * subsample_frac)
      # Ensure at least 2 per group
      n_sub_0 <- max(n_sub_0, 2)
      n_sub_1 <- max(n_sub_1, 2)
      boot_idx <- c(
        sample(idx_0, min(n_sub_0, length(idx_0)), replace = FALSE),
        sample(idx_1, min(n_sub_1, length(idx_1)), replace = FALSE)
      )
    } else {
      boot_idx <- sample(n, n_sub, replace = FALSE)
    }

    X_boot <- X[boot_idx, , drop = FALSE]
    y_boot <- y[boot_idx]

    # --- Run FS method -------------------------------------------------------
    result <- run_fs_method(method_name, X_boot, y_boot, params)

    selection_matrix[b, ]  <- as.integer(result$selected)
    importance_matrix[b, ] <- result$importance
    timings[b]    <- result$time
    converged[b]  <- result$converged
    messages[b]   <- result$message
  }

  list(
    method = method_name,
    category = method_config$category,
    selection_matrix = selection_matrix,
    importance_matrix = importance_matrix,
    timings = timings,
    converged = converged,
    messages = messages,
    n_boot = n_boot,
    n_converged = sum(converged),
    mean_n_selected = mean(rowSums(selection_matrix[converged, , drop = FALSE])),
    total_time = sum(timings)
  )
}

# ==============================================================================
# Main loop: iterate over datasets × methods
# ==============================================================================

cli::cli_h1("Running Feature Selection")

setup_parallel(config)
progressr::handlers(global = TRUE)

# Build task list (dataset × method), skipping completed checkpoints
tasks <- list()
for (sim_file in sim_files) {
  task_id <- gsub("^sim_|\\.rds$", "", basename(sim_file))

  for (method in methods) {
    fs_id <- paste0(task_id, "_", method$name)
    out_file <- file.path(fs_dir, paste0("fs_", fs_id, ".rds"))

    if (file.exists(out_file)) {
      next  # Already computed
    }

    tasks[[fs_id]] <- list(
      fs_id = fs_id,
      sim_file = sim_file,
      task_id = task_id,
      method = method
    )
  }
}

cli::cli_alert_info("FS tasks to run: {length(tasks)} (skipping existing)")

if (length(tasks) > 0) {
  # Process in chunks to manage memory
  chunk_size <- config$project$n_cores * 2
  task_list <- unname(tasks)
  n_chunks <- ceiling(length(task_list) / chunk_size)

  for (chunk_idx in seq_len(n_chunks)) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx   <- min(chunk_idx * chunk_size, length(task_list))
    chunk <- task_list[start_idx:end_idx]

    cli::cli_alert_info("Chunk {chunk_idx}/{n_chunks} ({length(chunk)} tasks)")

    # Parallelize at the dataset × method level
    chunk_results <- furrr::future_map(chunk, function(task) {
      # Load simulated data
      sim_data <- readRDS(task$sim_file)

      base_seed <- derive_seed(config$project$seed,
                               paste0("fs_", task$fs_id))

      # Run bootstrap FS
      result <- tryCatch(
        run_bootstrap_fs(
          sim_data = sim_data,
          method_config = task$method,
          n_boot = n_bootstrap,
          subsample_frac = subsample_frac,
          stratified = stratified,
          base_seed = base_seed
        ),
        error = function(e) {
          list(method = task$method$name, error = e$message,
               fs_id = task$fs_id)
        }
      )

      # Add metadata
      result$task_id <- task$task_id
      result$scenario <- sim_data$scenario
      result$level_label <- sim_data$level_label
      result$rep <- sim_data$rep
      result$true_features <- sim_data$true_features

      # Save checkpoint
      out_file <- file.path(fs_dir, paste0("fs_", task$fs_id, ".rds"))
      saveRDS(result, out_file)

      # Return summary
      list(fs_id = task$fs_id,
           method = task$method$name,
           status = if (is.null(result$error)) "success" else "error",
           n_converged = result$n_converged,
           mean_selected = round(result$mean_n_selected, 1),
           time = round(result$total_time, 1))
    }, .options = furrr::furrr_options(seed = TRUE),
       .progress = TRUE)

    # Log chunk results
    n_ok <- sum(vapply(chunk_results, function(r) r$status == "success", logical(1)))
    cli::cli_alert_success("  Chunk {chunk_idx}: {n_ok}/{length(chunk)} succeeded")

    # Periodic garbage collection
    gc(verbose = FALSE)
  }
}

teardown_parallel()

# ==============================================================================
# Summary
# ==============================================================================

cli::cli_h1("Feature Selection Summary")

fs_files <- list.files(fs_dir, pattern = "^fs_.*\\.rds$")
cli::cli_alert_info("Total FS result files: {length(fs_files)}")

# Summary per method
for (method in methods) {
  n_files <- sum(grepl(paste0("_", method$name, "\\.rds$"), fs_files))
  cli::cli_alert_info("  {method$name}: {n_files} results")
}

finalize_logging("Script 05", session_file = "results/session_info_05.txt")

cat("\n=== Script 05: Complete ===\n")
