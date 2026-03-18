# ==============================================================================
# Script 04 — Simulation framework
# ==============================================================================
# Purpose: Generate simulated metabolomics datasets for all scenarios defined
#          in config.yaml, using empirical parameters from real data.
#
# Design:  For each scenario × level × replication:
#          1. Generate correlated features from multivariate normal
#          2. Inject true signals (fold changes) in p_true features
#          3. Optionally add confounders and interactions
#          4. Introduce missing values (MCAR/MNAR)
#          5. Apply preprocessing
#          6. Save with ground truth labels
#
# Input:   data/empirical_params/empirical_params.rds
# Output:  data/simulated/sim_<scenario>_<level>_rep<r>.rds
# ==============================================================================

cat("=== Script 04: Simulation ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(mvtnorm)
  library(corpcor)
  library(MASS)
  library(future)
  library(furrr)
  library(progressr)
})

source(here("R", "utils", "helpers.R"))
config <- load_config()
setup_logging("Script 04: Simulation Framework", config)

params_dir <- get_path(config, "empirical_params")
sim_dir    <- get_path(config, "simulated_data")

# --- Load empirical parameters ------------------------------------------------
emp_params <- load_result(file.path(params_dir, "empirical_params.rds"))$data

# ==============================================================================
# Core simulation function
# ==============================================================================

#' Simulate a single metabolomics dataset
#'
#' @param n_case Number of case samples
#' @param n_control Number of control samples
#' @param p Total number of features
#' @param p_true Number of truly differential features
#' @param fc Fold change magnitude
#' @param fc_distribution How FC values are assigned ("fixed","uniform","decreasing")
#' @param cor_matrix Correlation matrix (p × p) or NULL
#' @param correlation_source Source of correlation structure ("empirical", "ar1", "block")
#' @param correlation_scale Scaling factor for off-diagonal correlations
#' @param missing_rate Proportion of missing values to introduce
#' @param missing_mechanism "MCAR" or "MNAR"
#' @param preprocessing Preprocessing method
#' @param n_confounders Number of confounding variables
#' @param n_interactions Number of feature interactions
#' @param seed Random seed
#' @return List: X, y, true_features, params
simulate_dataset <- function(n_case, n_control, p, p_true,
                             fc = 1.5, fc_distribution = "fixed",
                             cor_matrix = NULL, correlation_source = "empirical",
                             correlation_scale = 1.0,
                             missing_rate = 0.05, missing_mechanism = "MNAR",
                             preprocessing = "log_auto",
                             n_confounders = 0, n_interactions = 0,
                             seed = 42) {
  set.seed(seed)

  n <- n_case + n_control
  y <- c(rep(0L, n_control), rep(1L, n_case))

  # --- Step 1: Correlation matrix -------------------------------------------
  if (is.null(cor_matrix) || ncol(cor_matrix) != p) {
    if (identical(correlation_source, "ar1")) {
      rho <- 0.3 * correlation_scale
      cor_matrix <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
    } else if (identical(correlation_source, "block")) {
      block_size <- min(50L, p)
      rho <- 0.3 * correlation_scale
      cor_matrix <- diag(1, p)
      for (block_start in seq(1, p, by = block_size)) {
        block_end <- min(block_start + block_size - 1L, p)
        block_idx <- block_start:block_end
        block_cor <- matrix(rho, nrow = length(block_idx), ncol = length(block_idx))
        diag(block_cor) <- 1
        cor_matrix[block_idx, block_idx] <- block_cor
      }
    } else {
      cor_matrix <- build_correlation_matrix(p, emp_params, correlation_scale)
    }
  } else {
    # Scale off-diagonal elements
    if (correlation_scale != 1.0) {
      diag_vals <- diag(cor_matrix)
      cor_matrix <- cor_matrix * correlation_scale
      diag(cor_matrix) <- diag_vals
      # Ensure valid correlation matrix
      cor_matrix <- corpcor::make.positive.definite(cor_matrix)
      cor_matrix <- cov2cor(cor_matrix)
    }
  }

  # --- Step 2: Generate multivariate normal data ----------------------------
  mu <- rep(0, p)
  X <- mvtnorm::rmvnorm(n, mean = mu, sigma = cor_matrix)
  colnames(X) <- paste0("V", seq_len(p))

  # Transform to positive scale (simulate concentration data)
  # Exponentiate to get log-normal-like distributions
  X <- exp(X)

  # --- Step 3: Select true features and inject signal -----------------------
  true_idx <- sort(sample(seq_len(p), p_true))
  true_features <- rep(FALSE, p)
  true_features[true_idx] <- TRUE

  # Determine fold changes
  fc_values <- switch(fc_distribution,
    "fixed" = rep(fc, p_true),
    "uniform" = runif(p_true, min = 1.2, max = fc * 1.5),
    "decreasing" = seq(fc * 1.5, 1.2, length.out = p_true),
    rep(fc, p_true)
  )

  # Randomly assign direction (up/down regulation)
  directions <- sample(c(1, -1), p_true, replace = TRUE)

  # Apply fold changes to case samples
  for (i in seq_along(true_idx)) {
    j <- true_idx[i]
    if (directions[i] == 1) {
      X[y == 1, j] <- X[y == 1, j] * fc_values[i]
    } else {
      X[y == 1, j] <- X[y == 1, j] / fc_values[i]
    }
  }

  # --- Step 4: Add confounders ----------------------------------------------
  if (n_confounders > 0) {
    # Confounders: variables correlated with class but NOT true signal
    for (k in seq_len(n_confounders)) {
      # Confounder affects a random subset of non-true features
      conf_idx <- sample(which(!true_features), min(20, sum(!true_features)))
      conf_effect <- rnorm(1, mean = 0.3, sd = 0.1)
      # Confounder is partially correlated with class
      conf_var <- rbinom(n, 1, prob = ifelse(y == 1, 0.7, 0.3))
      for (j in conf_idx) {
        X[, j] <- X[, j] * (1 + conf_effect * conf_var)
      }
    }
  }

  # --- Step 5: Add interactions ---------------------------------------------
  if (n_interactions > 0) {
    # Create interaction effects between pairs of true features
    n_int_actual <- min(n_interactions, choose(p_true, 2))
    if (n_int_actual > 0 && p_true >= 2) {
      int_pairs <- combn(true_idx, 2)
      int_pairs <- int_pairs[, sample(ncol(int_pairs),
                                      min(n_int_actual, ncol(int_pairs))),
                             drop = FALSE]
      for (k in seq_len(ncol(int_pairs))) {
        j1 <- int_pairs[1, k]
        j2 <- int_pairs[2, k]
        # Interaction: effect depends on both features being high in cases
        int_effect <- X[, j1] * X[, j2] * 0.01
        X[y == 1, j1] <- X[y == 1, j1] + int_effect[y == 1]
      }
    }
  }

  # --- Step 6: Introduce missing values -------------------------------------
  if (missing_rate > 0) {
    n_missing <- round(n * p * missing_rate)
    if (missing_mechanism == "MCAR") {
      miss_idx <- sample(length(X), n_missing)
      X[miss_idx] <- NA
    } else if (missing_mechanism == "MNAR") {
      # MNAR: low-abundance values more likely to be missing
      # (common in LC-MS: below limit of detection)
      probs <- 1 / (1 + exp(scale(X)))  # Higher prob for lower values
      probs <- probs / sum(probs) * n_missing
      probs <- pmin(probs, 1)
      miss_mask <- matrix(rbinom(length(X), 1, prob = probs), nrow = n)
      X[miss_mask == 1] <- NA
    }
  }

  # --- Step 7: Impute missing values ----------------------------------------
  if (any(is.na(X))) {
    # Median imputation (matching Script 02)
    for (j in seq_len(ncol(X))) {
      na_idx <- is.na(X[, j])
      if (any(na_idx)) {
        X[na_idx, j] <- median(X[, j], na.rm = TRUE)
      }
    }
    X[is.na(X)] <- 0
  }

  # --- Step 8: Preprocess ----------------------------------------------------
  X_processed <- normalize_sim(X, method = preprocessing)

  # --- Return ----------------------------------------------------------------
  list(
    X = X_processed,
    X_raw = X,
    y = y,
    true_features = true_features,
    true_idx = true_idx,
    fc_values = fc_values,
    fc_directions = directions,
    params = list(
      n_case = n_case, n_control = n_control, p = p, p_true = p_true,
      fc = fc, fc_distribution = fc_distribution,
      correlation_source = correlation_source,
      correlation_scale = correlation_scale,
      missing_rate = missing_rate, missing_mechanism = missing_mechanism,
      preprocessing = preprocessing,
      n_confounders = n_confounders, n_interactions = n_interactions,
      seed = seed
    )
  )
}

#' Build a correlation matrix of dimension p from empirical data
#' @param p Target dimension
#' @param emp_params Empirical parameters list
#' @param scale Correlation scaling factor
#' @return p × p positive-definite correlation matrix
build_correlation_matrix <- function(p, emp_params, scale = 1.0) {
  # Use largest available empirical correlation matrix
  ref_cor <- NULL
  for (plat in names(emp_params$pooled)) {
    cm <- emp_params$pooled[[plat]]$cor_matrix
    if (!is.null(cm) && (is.null(ref_cor) || ncol(cm) > ncol(ref_cor))) {
      ref_cor <- cm
    }
  }

  if (!is.null(ref_cor) && ncol(ref_cor) >= p) {
    # Subsample from empirical
    idx <- sort(sample(ncol(ref_cor), p))
    sigma <- ref_cor[idx, idx]
  } else if (!is.null(ref_cor)) {
    # Tile / extend the empirical correlation
    p_emp <- ncol(ref_cor)
    sigma <- matrix(0, p, p)
    diag(sigma) <- 1
    # Fill with blocks from empirical
    for (i in seq(1, p, by = p_emp)) {
      j_end <- min(i + p_emp - 1, p)
      block_size <- j_end - i + 1
      sigma[i:j_end, i:j_end] <- ref_cor[1:block_size, 1:block_size]
    }
  } else {
    # Fallback: AR(1) structure
    cli::cli_alert_warning("No empirical correlation available, using AR(1)")
    rho <- 0.3 * scale
    sigma <- rho^abs(outer(1:p, 1:p, "-"))
  }

  # Apply scaling
  if (scale != 1.0 && !is.null(ref_cor)) {
    diag_vals <- diag(sigma)
    sigma <- sigma * scale
    diag(sigma) <- diag_vals
  }

  # Ensure positive definiteness
  sigma <- corpcor::make.positive.definite(sigma)
  sigma <- cov2cor(sigma)
  sigma
}

#' Normalize simulated data
#' @param X Matrix
#' @param method Method name
#' @return Processed matrix
normalize_sim <- function(X, method = "log_auto") {
  switch(method,
    "none" = X,
    "log_auto" = {
      X_log <- log2(pmax(X, 1e-10))
      scale(X_log, center = TRUE, scale = TRUE)
    },
    "log_pareto" = {
      X_log <- log2(pmax(X, 1e-10))
      X_centered <- scale(X_log, center = TRUE, scale = FALSE)
      sds <- apply(X_log, 2, sd, na.rm = TRUE)
      sweep(X_centered, 2, sqrt(pmax(sds, 1e-10)), "/")
    },
    "pqn_auto" = {
      ref <- apply(X, 2, median, na.rm = TRUE)
      quotients <- sweep(X, 2, pmax(ref, 1e-10), "/")
      dilution <- apply(quotients, 1, median, na.rm = TRUE)
      X_pqn <- sweep(X, 1, pmax(dilution, 1e-10), "/")
      X_log <- log2(pmax(X_pqn, 1e-10))
      scale(X_log, center = TRUE, scale = TRUE)
    },
    {
      X_log <- log2(pmax(X, 1e-10))
      scale(X_log, center = TRUE, scale = TRUE)
    }
  )
}

# ==============================================================================
# Generate all scenarios
# ==============================================================================

cli::cli_h1("Generating simulated datasets")

base <- config$simulation$base
scenarios <- config$simulation$scenarios
n_rep <- config$simulation$resampling$n_replications

setup_parallel(config)

# Build task list
tasks <- list()
for (scenario in scenarios) {
  for (level in scenario$levels) {
    for (rep in seq_len(n_rep)) {
      task_id <- paste0(scenario$name, "_", level$label, "_rep", rep)

      # Check for existing checkpoint
      if (checkpoint_exists(paste0("sim_", task_id), config)) {
        next
      }

      # Merge base params with level overrides
      sim_params <- base
      for (nm in names(level)) {
        if (nm != "label") sim_params[[nm]] <- level[[nm]]
      }

      tasks[[task_id]] <- list(
        task_id = task_id,
        scenario = scenario$name,
        level_label = level$label,
        rep = rep,
        params = sim_params
      )
    }
  }
}

cli::cli_alert_info("Tasks to simulate: {length(tasks)} (skipping existing checkpoints)")

if (length(tasks) > 0) {
  # Run simulations in parallel
  progressr::handlers(global = TRUE)

  results <- furrr::future_map(tasks, function(task) {
    seed <- derive_seed(config$project$seed,
                        paste0("sim_", task$task_id))
    sp <- task$params

    sim <- tryCatch(
      simulate_dataset(
        n_case = sp$n_case,
        n_control = sp$n_control,
        p = sp$p,
        p_true = sp$p_true,
        fc = sp$fc,
        fc_distribution = sp$fc_distribution,
        cor_matrix = NULL,
        correlation_source = sp$correlation_source,
        correlation_scale = sp$correlation_scale,
        missing_rate = sp$missing_rate,
        missing_mechanism = sp$missing_mechanism,
        preprocessing = sp$preprocessing,
        n_confounders = sp$n_confounders,
        n_interactions = sp$n_interactions,
        seed = seed
      ),
      error = function(e) {
        list(error = e$message, task_id = task$task_id)
      }
    )

    # Add task metadata
    sim$task_id <- task$task_id
    sim$scenario <- task$scenario
    sim$level_label <- task$level_label
    sim$rep <- task$rep

    # Save individual checkpoint
    fpath <- file.path(sim_dir, paste0("sim_", task$task_id, ".rds"))
    saveRDS(sim, fpath)

    # Return summary only (to avoid memory bloat)
    list(task_id = task$task_id, status = if (is.null(sim$error)) "success" else "error",
         n = length(sim$y), p = length(sim$true_features),
         p_true = sum(sim$true_features, na.rm = TRUE))
  }, .options = furrr::furrr_options(seed = TRUE),
     .progress = TRUE)

  cli::cli_alert_success("Simulation complete: {sum(vapply(results, function(r) r$status == 'success', logical(1)))} / {length(results)} succeeded")
} else {
  cli::cli_alert_info("All simulations already checkpointed")
}

teardown_parallel()

# ==============================================================================
# Summary
# ==============================================================================

cli::cli_h1("Simulation Summary")

# Count files per scenario
sim_files <- list.files(sim_dir, pattern = "^sim_.*\\.rds$")
cli::cli_alert_info("Total simulated datasets: {length(sim_files)}")

for (scenario in scenarios) {
  n_files <- sum(grepl(paste0("sim_", scenario$name), sim_files))
  n_expected <- length(scenario$levels) * n_rep
  cli::cli_alert_info("  {scenario$name}: {n_files}/{n_expected}")
}

finalize_logging("Script 04", session_file = "results/session_info_04.txt")

cat("\n=== Script 04: Complete ===\n")
