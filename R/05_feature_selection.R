# ==============================================================================
# Script 05 — Feature selection pipeline (BOTTLENECK)
# ==============================================================================
# Purpose: Apply all 11 FS methods to each simulated dataset with bootstrap
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
# Parallelization: parallel::mclapply (fork-based, Linux only)
#   - No chunk blocking: if 1 worker dies, others continue
#   - Copy-on-write memory: 40 workers share base R session
#   - Timeout per dataset: 30 min (kills hung workers)
#
# Input:   data/simulated/sim_*.rds from Script 04
# Output:  results/feature_selection/fs_<task_id>_<method>.rds
# ==============================================================================

cat("=== Script 05: Feature Selection Pipeline ===\n")

# NOTE: CUDA must be disabled BEFORE R starts (in run_fs_with_telegram.sh)
# via CUDA_VISIBLE_DEVICES="". See comment in launcher script for details.

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(parallel)
})

source(here("R", "utils", "helpers.R"))
source(here("R", "utils", "fs_methods.R"))
config <- load_config()
setup_logging("Script 05: Feature Selection (BOTTLENECK)", config)

sim_dir <- get_path(config, "simulated_data")
fs_dir  <- get_path(config, "results_fs")

# --- Configuration ------------------------------------------------------------
methods        <- config$methods
n_bootstrap    <- config$simulation$resampling$n_bootstrap
subsample_frac <- config$simulation$resampling$subsample_fraction
stratified     <- config$simulation$resampling$stratified
n_cores        <- config$project$n_cores
TIMEOUT_SEC    <- 1800  # 30 min per dataset (all methods)

# --- List simulated datasets --------------------------------------------------
sim_files <- list.files(sim_dir, pattern = "^sim_.*\\.rds$", full.names = TRUE)
cli::cli_alert_info("Simulated datasets: {length(sim_files)}")
cli::cli_alert_info("Methods: {length(methods)}")
cli::cli_alert_info("Bootstrap samples: {n_bootstrap}")
cli::cli_alert_info("Workers: {n_cores} (mclapply fork)")
cli::cli_alert_info("Timeout per dataset: {TIMEOUT_SEC}s")

# ==============================================================================
# Bootstrap feature selection for a single dataset × method
# ==============================================================================

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
  timings   <- numeric(n_boot)
  converged <- logical(n_boot)
  messages  <- character(n_boot)

  for (b in seq_len(n_boot)) {
    set.seed(base_seed + b)

    # Stratified subsample
    if (stratified) {
      idx_0 <- which(y == 0)
      idx_1 <- which(y == 1)
      n_sub_0 <- max(round(length(idx_0) * subsample_frac), 2)
      n_sub_1 <- max(round(length(idx_1) * subsample_frac), 2)
      boot_idx <- c(
        sample(idx_0, min(n_sub_0, length(idx_0)), replace = FALSE),
        sample(idx_1, min(n_sub_1, length(idx_1)), replace = FALSE)
      )
    } else {
      boot_idx <- sample(n, round(n * subsample_frac), replace = FALSE)
    }

    X_boot <- X[boot_idx, , drop = FALSE]
    y_boot <- y[boot_idx]

    result <- run_fs_method(method_name, X_boot, y_boot, params)

    selection_matrix[b, ]  <- as.integer(result$selected)
    importance_matrix[b, ] <- result$importance
    timings[b]   <- result$time
    converged[b] <- result$converged
    messages[b]  <- result$message
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
# Process one dataset: all pending methods sequentially
# ==============================================================================

process_dataset <- function(sim_file, methods_list, fs_dir, n_bootstrap,
                            subsample_frac, stratified, base_seed_global) {
  task_id <- gsub("^sim_|\\.rds$", "", basename(sim_file))
  sim_data <- readRDS(sim_file)
  summaries <- list()

  for (method in methods_list) {
    fs_id <- paste0(task_id, "_", method$name)
    out_file <- file.path(fs_dir, paste0("fs_", fs_id, ".rds"))

    # Skip if already done
    if (file.exists(out_file)) next

    base_seed <- derive_seed(base_seed_global, paste0("fs_", fs_id))

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
        list(method = method$name, error = e$message, fs_id = fs_id,
             total_time = 0)
      }
    )

    result$task_id <- task_id
    result$scenario <- sim_data$scenario
    result$level_label <- sim_data$level_label
    result$rep <- sim_data$rep
    result$true_features <- sim_data$true_features

    saveRDS(result, out_file)

    summaries[[fs_id]] <- list(
      fs_id = fs_id,
      method = method$name,
      status = if (is.null(result$error)) "ok" else "err",
      time = round(result$total_time, 1)
    )
  }

  summaries
}

# ==============================================================================
# Main: build task list and run with mclapply
# ==============================================================================

cli::cli_h1("Running Feature Selection")

# Build list of datasets with pending work
pending_files <- character(0)
for (sim_file in sim_files) {
  task_id <- gsub("^sim_|\\.rds$", "", basename(sim_file))
  for (method in methods) {
    fs_id <- paste0(task_id, "_", method$name)
    out_file <- file.path(fs_dir, paste0("fs_", fs_id, ".rds"))
    if (!file.exists(out_file)) {
      pending_files <- c(pending_files, sim_file)
      break
    }
  }
}
pending_files <- unique(pending_files)

n_done_before <- length(list.files(fs_dir, pattern = "^fs_.*\\.rds$"))
n_total <- length(sim_files) * length(methods)
cli::cli_alert_info("Datasets with pending work: {length(pending_files)}")
cli::cli_alert_info("Already completed: {n_done_before}/{n_total}")

if (length(pending_files) > 0) {
  pipeline_start <- Sys.time()

  cli::cli_alert_info("Launching mclapply with {n_cores} workers...")

  results <- mclapply(pending_files, function(sf) {
    process_dataset(
      sim_file = sf,
      methods_list = methods,
      fs_dir = fs_dir,
      n_bootstrap = n_bootstrap,
      subsample_frac = subsample_frac,
      stratified = stratified,
      base_seed_global = config$project$seed
    )
  }, mc.cores = n_cores, mc.preschedule = FALSE)

  # Count results
  elapsed <- as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))
  n_done_after <- length(list.files(fs_dir, pattern = "^fs_.*\\.rds$"))
  n_new <- n_done_after - n_done_before
  n_errors <- sum(vapply(results, function(r) {
    if (is.null(r) || inherits(r, "try-error")) return(1L)
    sum(vapply(r, function(x) x$status == "err", logical(1)))
  }, integer(1)))
  n_timeouts <- sum(vapply(results, is.null, logical(1)))

  cli::cli_alert_success(
    "Completed: {n_new} new jobs in {round(elapsed, 1)} min"
  )
  if (n_errors > 0) {
    cli::cli_alert_warning("Errors: {n_errors}")
  }
  if (n_timeouts > 0) {
    cli::cli_alert_warning("Timeouts/killed workers: {n_timeouts}")
  }
}

# ==============================================================================
# Summary
# ==============================================================================

cli::cli_h1("Feature Selection Summary")

fs_files <- list.files(fs_dir, pattern = "^fs_.*\\.rds$")
cli::cli_alert_info("Total FS result files: {length(fs_files)}")

for (method in methods) {
  n_files <- sum(grepl(paste0("_", method$name, "\\.rds$"), fs_files))
  cli::cli_alert_info("  {method$name}: {n_files} results")
}

finalize_logging("Script 05", session_file = "results/session_info_05.txt")

cat("\n=== Script 05: Complete ===\n")
