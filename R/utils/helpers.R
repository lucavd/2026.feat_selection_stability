# ==============================================================================
# Utility functions — helpers.R
# ==============================================================================
# Shared functions used across all scripts.
# Source with: source(here::here("R", "utils", "helpers.R"))
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(cli)
  library(fs)
  library(data.table)
})

# --- Configuration loader -----------------------------------------------------

#' Load project configuration
#' @return Named list from config.yaml
load_config <- function() {
  cfg_path <- here::here("config", "config.yaml")
  if (!file.exists(cfg_path)) {
    cli::cli_abort("Configuration file not found: {cfg_path}")
  }
  yaml::read_yaml(cfg_path)
}

# --- Reproducibility ----------------------------------------------------------

#' Derive a deterministic seed from base seed and a string identifier
#' @param base_seed Integer base seed from config
#' @param id Character string (e.g. scenario name, replicate number)
#' @return Integer seed
derive_seed <- function(base_seed, id) {
  hash <- digest::digest(paste0(base_seed, "_", id), algo = "xxhash32",
                         serialize = FALSE)
  # Convert hex hash to integer, mod 2^31 for R compatibility
  as.integer(abs(strtoi(substr(hash, 1, 7), base = 16L)) %% .Machine$integer.max)
}

#' Set seed reproducibly for a given operation
#' @param config Config list (must contain project$seed)
#' @param id Character identifier for this operation
set_project_seed <- function(config, id) {
  seed <- derive_seed(config$project$seed, id)
  set.seed(seed)
  invisible(seed)
}

# --- Path helpers -------------------------------------------------------------

#' Resolve a path from config, creating directory if needed
#' @param config Config list
#' @param key Character key in config$paths (e.g. "raw_data")
#' @return Absolute path
get_path <- function(config, key) {
  rel <- config$paths[[key]]
  if (is.null(rel)) {
    cli::cli_abort("Path key '{key}' not found in config$paths")
  }
  p <- here::here(rel)
  fs::dir_create(p)
  p
}

# --- Logging ------------------------------------------------------------------

#' Initialize logging for a script
#' @param script_name Character name of script
#' @param config Config list
setup_logging <- function(script_name, config = NULL) {
  cli::cli_h1("{script_name}")
  cli::cli_alert_info("Start time: {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}")
  if (!is.null(config)) {
    cli::cli_alert_info("Seed: {config$project$seed}")
    cli::cli_alert_info("Cores: {config$project$n_cores}")
  }
}

#' Finalize logging for a script
#' @param script_name Character name of script
#' @param session_file Optional file to save sessionInfo
finalize_logging <- function(script_name, session_file = NULL) {
  cli::cli_alert_success("{script_name} completed at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}")
  if (!is.null(session_file)) {
    writeLines(capture.output(sessionInfo()), here::here(session_file))
    cli::cli_alert_info("Session info saved to {session_file}")
  }
}

# --- Checkpointing -----------------------------------------------------------

#' Save a checkpoint (RDS)
#' @param data Object to save
#' @param name Character checkpoint name
#' @param config Config list
save_checkpoint <- function(data, name, config) {
  cp_dir <- get_path(config, "checkpoints")
  fpath  <- file.path(cp_dir, paste0(name, ".rds"))
  saveRDS(data, fpath)
  cli::cli_alert_success("Checkpoint saved: {name}")
  invisible(fpath)
}

#' Load a checkpoint if it exists
#' @param name Character checkpoint name
#' @param config Config list
#' @return Object or NULL
load_checkpoint <- function(name, config) {
  cp_dir <- get_path(config, "checkpoints")
  fpath  <- file.path(cp_dir, paste0(name, ".rds"))
  if (file.exists(fpath)) {
    cli::cli_alert_info("Resuming from checkpoint: {name}")
    return(readRDS(fpath))
  }
  NULL
}

#' Check if a checkpoint exists
#' @param name Character checkpoint name
#' @param config Config list
#' @return Logical
checkpoint_exists <- function(name, config) {
  cp_dir <- here::here(config$paths$checkpoints)
  file.exists(file.path(cp_dir, paste0(name, ".rds")))
}

# --- Data I/O -----------------------------------------------------------------

#' Save result as RDS with metadata
#' @param data Object to save
#' @param path File path
#' @param metadata Named list of metadata
save_result <- function(data, path, metadata = list()) {
  result <- list(
    data = data,
    metadata = c(metadata, list(
      timestamp = Sys.time(),
      r_version = R.version.string
    ))
  )
  saveRDS(result, path)
  cli::cli_alert_success("Saved: {basename(path)}")
  invisible(path)
}

#' Load result RDS
#' @param path File path
#' @return List with data and metadata
load_result <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("File not found: {path}")
  }
  readRDS(path)
}

# --- Parallelization setup ---------------------------------------------------

#' Setup future plan for parallel execution
#' @param config Config list
setup_parallel <- function(config) {
  n_cores_requested <- config$project$n_cores
  n_cores <- min(n_cores_requested, 100L)
  future::plan(future::multisession, workers = n_cores)
  if (n_cores < n_cores_requested) {
    cli::cli_alert_warning("Requested {n_cores_requested} workers; capping at {n_cores}")
  }
  cli::cli_alert_info("Parallel plan: multisession with {n_cores} workers")
}

#' Teardown parallel plan
teardown_parallel <- function() {
  future::plan(future::sequential)
  cli::cli_alert_info("Parallel plan reset to sequential")
}

# --- Validation ---------------------------------------------------------------

#' Validate a feature matrix
#' @param X Matrix or data.frame (samples × features)
#' @param y Response vector
#' @param label Character label for messages
validate_data <- function(X, y, label = "data") {
  stopifnot(is.matrix(X) || is.data.frame(X))
  stopifnot(length(y) == nrow(X))
  n <- nrow(X)
  p <- ncol(X)
  n_na <- sum(is.na(X))
  pct_na <- round(100 * n_na / (n * p), 1)
  cli::cli_alert_info("{label}: n={n}, p={p}, missing={pct_na}%")
  invisible(TRUE)
}
