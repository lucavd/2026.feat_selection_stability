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

# Build dataset-level task list: one task per dataset, all methods inside
dataset_tasks <- list()
for (sim_file in sim_files) {
  task_id <- gsub("^sim_|\\.rds$", "", basename(sim_file))

  # Check which methods still need to run for this dataset
  pending_methods <- list()
  for (method in methods) {
    fs_id <- paste0(task_id, "_", method$name)
    out_file <- file.path(fs_dir, paste0("fs_", fs_id, ".rds"))
    if (!file.exists(out_file)) {
      pending_methods <- c(pending_methods, list(method))
    }
  }

  if (length(pending_methods) > 0) {
    dataset_tasks[[task_id]] <- list(
      task_id = task_id,
      sim_file = sim_file,
      methods = pending_methods
    )
  }
}

n_pending_jobs <- sum(vapply(dataset_tasks, function(dt) length(dt$methods), integer(1)))
cli::cli_alert_info("Datasets with pending work: {length(dataset_tasks)}")
cli::cli_alert_info("Total method×dataset jobs: {n_pending_jobs} (skipping existing)")

if (length(dataset_tasks) > 0) {
  # Process in chunks of n_cores datasets at a time
  chunk_size <- config$project$n_cores
  task_list <- unname(dataset_tasks)
  n_chunks <- ceiling(length(task_list) / chunk_size)

  for (chunk_idx in seq_len(n_chunks)) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx   <- min(chunk_idx * chunk_size, length(task_list))
    chunk <- task_list[start_idx:end_idx]

    n_jobs_chunk <- sum(vapply(chunk, function(dt) length(dt$methods), integer(1)))
    cli::cli_alert_info("Chunk {chunk_idx}/{n_chunks} ({length(chunk)} datasets, {n_jobs_chunk} jobs)")

    # Each worker processes one dataset: all methods sequentially
    chunk_results <- furrr::future_map(chunk, function(dt) {
      sim_data <- readRDS(dt$sim_file)
      results_summary <- list()

      for (method in dt$methods) {
        fs_id <- paste0(dt$task_id, "_", method$name)

        base_seed <- derive_seed(config$project$seed,
                                 paste0("fs_", fs_id))

        result <- tryCatch(
          run_bootstrap_fs(
            sim_data = sim_data,
            method_config = method,
            n_boot = n_bootstrap,
            subsample_frac = subsample_frac,
            stratified = stratified,
            base_seed = base_seed
          ),
          error = function(e) {
            list(method = method$name, error = e$message, fs_id = fs_id)
          }
        )

        result$task_id <- dt$task_id
        result$scenario <- sim_data$scenario
        result$level_label <- sim_data$level_label
        result$rep <- sim_data$rep
        result$true_features <- sim_data$true_features

        out_file <- file.path(fs_dir, paste0("fs_", fs_id, ".rds"))
        saveRDS(result, out_file)

        results_summary[[fs_id]] <- list(
          fs_id = fs_id, method = method$name,
          status = if (is.null(result$error)) "success" else "error",
          time = round(result$total_time, 1))
      }

      results_summary
    }, .options = furrr::furrr_options(seed = TRUE),
       .progress = TRUE)

    # Flatten and count
    all_summaries <- do.call(c, chunk_results)
    n_ok <- sum(vapply(all_summaries, function(r) r$status == "success", logical(1)))
    cli::cli_alert_success("  Chunk {chunk_idx}: {n_ok}/{length(all_summaries)} succeeded")

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
