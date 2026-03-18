# ==============================================================================
# Script 02 — Preprocess and harmonize datasets
# ==============================================================================
# Purpose: Parse raw downloads, apply QC filters, impute, normalize,
#          and produce harmonized feature matrices for each dataset.
# Input:   data/raw/<accession>/ from Script 01
# Output:  data/processed/<accession>_processed.rds
#          data/processed/datasets_summary.rds
# ==============================================================================

cat("=== Script 02: Preprocessing ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

source(here("R", "utils", "helpers.R"))
config <- load_config()
setup_logging("Script 02: Preprocessing", config)

raw_dir  <- get_path(config, "raw_data")
proc_dir <- get_path(config, "processed_data")

# Load download manifest
manifest <- load_result(file.path(raw_dir, "download_manifest.rds"))$data
successful <- names(manifest)[vapply(manifest, function(m) m$status == "success",
                                     logical(1))]

cli::cli_alert_info("{length(successful)} datasets available for preprocessing")

get_matching_colnames <- function(df, pattern) {
  idx <- grep(pattern, tolower(names(df)))
  if (length(idx) == 0) {
    character(0)
  } else {
    names(df)[idx]
  }
}

align_sample_info <- function(sample_info, sample_ids) {
  if (is.null(sample_info) || length(sample_ids) == 0) {
    return(list(sample_info = sample_info,
                matched = rep(NA_integer_, length(sample_ids))))
  }

  id_cols <- get_matching_colnames(
    sample_info,
    "sample[ _]?name|sample[ _]?id|subject[ _]?id|specimen[ _]?id|participant[ _]?id|sample"
  )

  if (length(id_cols) > 0) {
    matched <- match(sample_ids, as.character(sample_info[[id_cols[1]]]))
    if (sum(!is.na(matched)) > 0) {
      return(list(sample_info = sample_info[matched, , drop = FALSE],
                  matched = matched))
    }
  }

  cli::cli_alert_warning("Unable to align sample metadata by ID; class labels will be set to NA for this dataset")

  list(sample_info = sample_info,
       matched = rep(NA_integer_, length(sample_ids)))
}

# ==============================================================================
# Parsing functions
# ==============================================================================

#' Parse MetaboLights MAF (Metabolite Assignment File) + sample sheet
#' @param acc_id Accession ID
#' @param raw_dir Raw data directory
#' @return List with X (matrix), y (vector), sample_info, feature_info
parse_metabolights <- function(acc_id, raw_dir) {
  study_dir <- file.path(raw_dir, acc_id)

  # Find MAF file (metabolite data)
  maf_files <- list.files(study_dir, pattern = "^m_.*\\.tsv$", full.names = TRUE)
  if (length(maf_files) == 0) {
    cli::cli_alert_danger("{acc_id}: no MAF file found")
    return(NULL)
  }

  # Find sample sheet
  sample_files <- list.files(study_dir, pattern = "^s_.*\\.txt$", full.names = TRUE)

  # Read MAF
  maf <- data.table::fread(maf_files[1], header = TRUE, check.names = FALSE)
  cli::cli_alert_info("MAF: {nrow(maf)} rows × {ncol(maf)} cols")

  # Identify metadata vs data columns
  # MetaboLights MAF: first columns are metabolite annotations,
  # then sample columns with numeric values
  meta_cols <- c("database_identifier", "chemical_formula", "smiles", "inchi",
                 "metabolite_identification", "mass_to_charge", "fragmentation",
                 "modifications", "charge", "retention_time",
                 "taxid", "species", "database", "database_version",
                 "reliability", "uri", "search_engine", "search_engine_score",
                 "smallmolecule_abundance_sub", "smallmolecule_abundance_stdev_sub",
                 "smallmolecule_abundance_std_error_sub")

  # Sample columns = all columns that are NOT metadata
  all_cols <- names(maf)
  potential_data_cols <- setdiff(all_cols, all_cols[tolower(all_cols) %in% tolower(meta_cols)])

  # Try to convert to numeric to identify actual data columns
  data_cols <- character(0)
  for (col in potential_data_cols) {
    vals <- suppressWarnings(as.numeric(maf[[col]]))
    if (sum(!is.na(vals)) > nrow(maf) * 0.5) {
      data_cols <- c(data_cols, col)
    }
  }

  if (length(data_cols) < 5) {
    cli::cli_alert_danger("{acc_id}: too few numeric columns ({length(data_cols)})")
    return(NULL)
  }

  # Build feature matrix (samples × features)
  X <- t(as.matrix(maf[, data_cols, with = FALSE]))
  X <- apply(X, 2, as.numeric)

  # Feature names from metabolite_identification or row index
  if ("metabolite_identification" %in% names(maf)) {
    feat_names <- make.unique(as.character(maf$metabolite_identification))
  } else {
    feat_names <- paste0("F", seq_len(ncol(X)))
  }
  colnames(X) <- feat_names
  rownames(X) <- data_cols

  # Parse sample sheet for class labels
  y <- NULL
  sample_info <- NULL
  if (length(sample_files) > 0) {
    samp <- data.table::fread(sample_files[1], header = TRUE, check.names = FALSE)
    factor_cols <- get_matching_colnames(
      samp,
      "factor|group|class|disease|condition|diagnosis|status"
    )
    if (length(factor_cols) > 0) {
      aligned <- align_sample_info(samp, rownames(X))
      sample_info <- aligned$sample_info
      if (sum(!is.na(aligned$matched)) > 0) {
        y_raw <- sample_info[[factor_cols[1]]]
        y <- as.character(y_raw)
      }
    }
  }

  feature_info <- if ("metabolite_identification" %in% names(maf)) {
    data.frame(
      feature_id = feat_names,
      original_name = maf$metabolite_identification,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(feature_id = feat_names, stringsAsFactors = FALSE)
  }

  list(X = X, y = y, sample_info = sample_info, feature_info = feature_info)
}

#' Parse Metabolomics Workbench JSON data
#' @param acc_id Accession ID
#' @param raw_dir Raw data directory
#' @return List with X (matrix), y (vector), sample_info, feature_info
parse_mw <- function(acc_id, raw_dir) {
  study_dir <- file.path(raw_dir, acc_id)

  # Try data matrix
  data_file <- file.path(study_dir, paste0(acc_id, "_data.json"))
  factors_file <- file.path(study_dir, paste0(acc_id, "_factors.json"))

  if (!file.exists(data_file)) {
    cli::cli_alert_danger("{acc_id}: no data file found")
    return(NULL)
  }

  # Parse data JSON
  data_json <- jsonlite::fromJSON(data_file, simplifyDataFrame = TRUE)

  # MW returns data in various formats; handle the common ones
  X <- NULL
  if (is.data.frame(data_json)) {
    # Direct data frame
    num_cols <- vapply(data_json, is.numeric, logical(1))
    if (sum(num_cols) > 5) {
      X <- as.matrix(data_json[, num_cols])
      sample_ids <- if ("sample_id" %in% names(data_json)) {
        data_json$sample_id
      } else {
        paste0("S", seq_len(nrow(X)))
      }
      rownames(X) <- sample_ids
    }
  } else if (is.list(data_json)) {
    # Nested format — try to extract
    cli::cli_alert_info("Nested JSON format, attempting extraction...")
    # Common MW pattern: list of metabolite -> sample -> value
    if (length(data_json) > 0) {
      mat <- tryCatch({
        df <- as.data.frame(data_json, check.names = FALSE)
        num_cols <- vapply(df, function(x) is.numeric(x) || all(!is.na(suppressWarnings(as.numeric(x)))),
                           logical(1))
        as.matrix(apply(df[, num_cols, drop = FALSE], 2, as.numeric))
      }, error = function(e) NULL)
      if (!is.null(mat) && ncol(mat) > 5) {
        X <- mat
      }
    }
  }

  if (is.null(X)) {
    cli::cli_alert_danger("{acc_id}: could not parse feature matrix")
    return(NULL)
  }

  # Parse factors for class labels
  y <- NULL
  sample_info <- NULL
  if (file.exists(factors_file)) {
    factors_json <- tryCatch(
      jsonlite::fromJSON(factors_file, simplifyDataFrame = TRUE),
      error = function(e) NULL
    )
    if (!is.null(factors_json) && is.data.frame(factors_json)) {
      class_cols <- get_matching_colnames(
        factors_json,
        "factor|group|class|disease|condition|diagnosis|status"
      )
      aligned <- align_sample_info(factors_json, rownames(X))
      sample_info <- aligned$sample_info
      if (length(class_cols) > 0 && sum(!is.na(aligned$matched)) > 0) {
        y <- as.character(sample_info[[class_cols[1]]])
      }
    }
  }

  feature_info <- data.frame(
    feature_id = colnames(X),
    stringsAsFactors = FALSE
  )

  list(X = X, y = y, sample_info = sample_info, feature_info = feature_info)
}

# ==============================================================================
# Quality Control and Filtering
# ==============================================================================

#' Apply QC filters to a dataset
#' @param dat List with X, y from parsing
#' @param config Config list
#' @return Filtered list or NULL if fails QC
apply_qc_filters <- function(dat, config) {
  X <- dat$X
  y <- dat$y
  ac <- config$acceptance_criteria

  n_orig <- nrow(X)
  p_orig <- ncol(X)

  # Step 1: Remove features with too many missing values
  feat_miss <- colMeans(is.na(X))
  keep_feat <- feat_miss <= ac$max_missing_rate_feature
  X <- X[, keep_feat, drop = FALSE]
  if (!is.null(dat$feature_info) && nrow(dat$feature_info) == length(keep_feat)) {
    dat$feature_info <- dat$feature_info[keep_feat, , drop = FALSE]
  }
  cli::cli_alert_info("Features: {p_orig} → {ncol(X)} (removed {sum(!keep_feat)} with >{ac$max_missing_rate_feature*100}% missing)")

  # Step 2: Remove samples with too many missing values
  samp_miss <- rowMeans(is.na(X))
  keep_samp <- samp_miss <= ac$max_missing_rate_sample
  X <- X[keep_samp, , drop = FALSE]
  if (!is.null(y)) y <- y[keep_samp]
  if (!is.null(dat$sample_info) && nrow(dat$sample_info) == length(keep_samp)) {
    dat$sample_info <- dat$sample_info[keep_samp, , drop = FALSE]
  }
  cli::cli_alert_info("Samples: {n_orig} → {nrow(X)} (removed {sum(!keep_samp)} with >{ac$max_missing_rate_sample*100}% missing)")

  # Step 3: Check minimum features
  if (ncol(X) < ac$min_features) {
    cli::cli_alert_danger("Too few features: {ncol(X)} < {ac$min_features}")
    return(NULL)
  }

  # Step 4: Remove zero-variance features
  feat_var <- apply(X, 2, var, na.rm = TRUE)
  keep_var <- !is.na(feat_var) & feat_var > 1e-10
  X <- X[, keep_var, drop = FALSE]
  if (!is.null(dat$feature_info) && nrow(dat$feature_info) == length(keep_var)) {
    dat$feature_info <- dat$feature_info[keep_var, , drop = FALSE]
  }
  cli::cli_alert_info("After zero-variance filter: {ncol(X)} features")

  dat$X <- X
  dat$y <- y
  dat
}

#' Binarize class labels
#' @param y Character vector of class labels
#' @return Named list: y_binary (0/1 integer), class_map, valid
binarize_classes <- function(y) {
  if (is.null(y)) return(list(y_binary = NULL, class_map = NULL, valid = FALSE))

  # Remove NA
  tab <- table(y[!is.na(y)])
  tab <- sort(tab, decreasing = TRUE)

  if (length(tab) < 2) {
    cli::cli_alert_warning("Only one class found: {names(tab)}")
    return(list(y_binary = NULL, class_map = NULL, valid = FALSE))
  }

  # Take top 2 classes
  top2 <- names(tab)[1:2]
  cli::cli_alert_info("Classes: '{top2[1]}' (n={tab[top2[1]]}) vs '{top2[2]}' (n={tab[top2[2]]})")

  # Assign: larger group = 0 (control), smaller = 1 (case)
  y_binary <- ifelse(y == top2[2], 1L, ifelse(y == top2[1], 0L, NA_integer_))

  list(
    y_binary = y_binary,
    class_map = c(control = top2[1], case = top2[2]),
    valid = TRUE
  )
}

#' Impute missing values (KNN or median)
#' @param X Numeric matrix
#' @return Imputed matrix
impute_missing <- function(X) {
  n_missing <- sum(is.na(X))
  if (n_missing == 0) {
    cli::cli_alert_info("No missing values to impute")
    return(X)
  }

  pct_missing <- round(100 * n_missing / length(X), 1)
  cli::cli_alert_info("Imputing {n_missing} values ({pct_missing}%)")

  # Use column median imputation (fast, conservative)
  for (j in seq_len(ncol(X))) {
    na_idx <- is.na(X[, j])
    if (any(na_idx)) {
      X[na_idx, j] <- median(X[, j], na.rm = TRUE)
    }
  }

  # Replace any remaining NAs with 0 (shouldn't happen)
  X[is.na(X)] <- 0
  X
}

#' Normalize a feature matrix
#' @param X Numeric matrix (samples × features)
#' @param method Normalization method from config
#' @return Normalized matrix
normalize_matrix <- function(X, method = "log_auto") {
  cli::cli_alert_info("Normalization: {method}")

  switch(method,
    "none" = X,

    "log_auto" = {
      # Log2 transform + autoscaling (mean-center, unit variance)
      X_log <- log2(pmax(X, 1))  # Avoid log(0)
      scale(X_log, center = TRUE, scale = TRUE)
    },

    "log_pareto" = {
      # Log2 transform + Pareto scaling (divide by sqrt(sd))
      X_log <- log2(pmax(X, 1))
      X_centered <- scale(X_log, center = TRUE, scale = FALSE)
      sds <- apply(X_log, 2, sd, na.rm = TRUE)
      sweep(X_centered, 2, sqrt(pmax(sds, 1e-10)), "/")
    },

    "pqn_auto" = {
      # Probabilistic Quotient Normalization + autoscaling
      ref <- apply(X, 2, median, na.rm = TRUE)  # Reference spectrum
      quotients <- sweep(X, 2, pmax(ref, 1e-10), "/")
      dilution <- apply(quotients, 1, median, na.rm = TRUE)
      X_pqn <- sweep(X, 1, pmax(dilution, 1e-10), "/")
      X_log <- log2(pmax(X_pqn, 1))
      scale(X_log, center = TRUE, scale = TRUE)
    },

    # Default fallback
    {
      cli::cli_alert_warning("Unknown method '{method}', using log_auto")
      X_log <- log2(pmax(X, 1))
      scale(X_log, center = TRUE, scale = TRUE)
    }
  )
}

# ==============================================================================
# Main preprocessing loop
# ==============================================================================

cli::cli_h1("Processing datasets")

datasets_summary <- list()

for (acc_id in successful) {
  cli::cli_h2("Processing {acc_id}")

  info <- manifest[[acc_id]]

  # Parse based on source
  parsed <- tryCatch({
    if (info$source == "MetaboLights") {
      parse_metabolights(acc_id, raw_dir)
    } else {
      parse_mw(acc_id, raw_dir)
    }
  }, error = function(e) {
    cli::cli_alert_danger("Parse error: {e$message}")
    NULL
  })

  if (is.null(parsed)) {
    cli::cli_alert_warning("{acc_id}: parsing failed, skipping")
    datasets_summary[[acc_id]] <- list(id = acc_id, status = "parse_failed")
    next
  }

  # Apply QC
  filtered <- apply_qc_filters(parsed, config)
  if (is.null(filtered)) {
    cli::cli_alert_warning("{acc_id}: failed QC, skipping")
    datasets_summary[[acc_id]] <- list(id = acc_id, status = "qc_failed")
    next
  }

  # Binarize classes
  class_result <- binarize_classes(filtered$y)
  if (!class_result$valid) {
    cli::cli_alert_warning("{acc_id}: cannot binarize classes, skipping")
    datasets_summary[[acc_id]] <- list(id = acc_id, status = "class_failed")
    next
  }

  # Keep only samples with valid class labels
  valid_idx <- !is.na(class_result$y_binary)
  X <- filtered$X[valid_idx, , drop = FALSE]
  y <- class_result$y_binary[valid_idx]
  sample_info <- filtered$sample_info
  if (!is.null(sample_info) && nrow(sample_info) == length(valid_idx)) {
    sample_info <- sample_info[valid_idx, , drop = FALSE]
  }

  # Check minimum samples per group
  tab <- table(y)
  if (any(tab < config$acceptance_criteria$min_samples_per_group)) {
    cli::cli_alert_warning("{acc_id}: insufficient samples per group ({paste(tab, collapse='/')}), skipping")
    datasets_summary[[acc_id]] <- list(id = acc_id, status = "too_few_samples")
    next
  }

  # Impute
  X_pre_impute <- X
  X <- impute_missing(X)

  # Normalize (using default preprocessing for now; Script 04 will vary this)
  X_norm <- normalize_matrix(X, method = config$simulation$base$preprocessing)

  # Final validation
  validate_data(X_norm, y, label = acc_id)

  # Save
  processed <- list(
    X = X_norm,
    X_raw = X_pre_impute,
    X_imputed = X,
    y = y,
    class_map = class_result$class_map,
    feature_info = filtered$feature_info,
    sample_info = sample_info,
    source = info$source,
    platform = info$platform,
    accession = acc_id,
    description = info$description
  )

  out_file <- file.path(proc_dir, paste0(acc_id, "_processed.rds"))
  save_result(processed, out_file,
              metadata = list(n = nrow(X_norm), p = ncol(X_norm),
                              n_case = sum(y == 1), n_control = sum(y == 0)))

  datasets_summary[[acc_id]] <- list(
    id = acc_id, status = "success",
    source = info$source, platform = info$platform,
    n = nrow(X_norm), p = ncol(X_norm),
    n_case = sum(y == 1), n_control = sum(y == 0)
  )

  cli::cli_alert_success("{acc_id}: {nrow(X_norm)} samples × {ncol(X_norm)} features")
}

# ==============================================================================
# Summary
# ==============================================================================

cli::cli_h1("Preprocessing Summary")

summary_df <- do.call(rbind, lapply(datasets_summary, function(s) {
  data.frame(
    id = s$id,
    status = s$status,
    source = ifelse(is.null(s$source), NA, s$source),
    platform = ifelse(is.null(s$platform), NA, s$platform),
    n = ifelse(is.null(s$n), NA, s$n),
    p = ifelse(is.null(s$p), NA, s$p),
    stringsAsFactors = FALSE
  )
}))

print(summary_df)

save_result(datasets_summary,
            file.path(proc_dir, "datasets_summary.rds"),
            metadata = list(script = "02_preprocess.R",
                            n_processed = sum(summary_df$status == "success", na.rm = TRUE)))

finalize_logging("Script 02", session_file = "results/session_info_02.txt")

cat("\n=== Script 02: Complete ===\n")
