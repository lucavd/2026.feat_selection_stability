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
  library(readxl)
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

#' Parse CIMCB benchmark Excel files (Mendez et al. 2019)
#' Format: Sheet "Data" with SampleID, Class, M1..Mn
#'         Sheet "Peak" with metabolite names
#' @param acc_id Accession ID
#' @param raw_dir Raw data directory
#' @param acc_config Config entry for this accession
#' @return List with X, y, sample_info, feature_info
parse_cimcb <- function(acc_id, raw_dir, acc_config) {
  xlsx_file <- file.path(raw_dir, acc_id, paste0(acc_id, ".xlsx"))
  if (!file.exists(xlsx_file)) {
    cli::cli_alert_danger("{acc_id}: Excel file not found")
    return(NULL)
  }

  data_df <- readxl::read_excel(xlsx_file, sheet = "Data")
  cli::cli_alert_info(
    "CIMCB Excel: {nrow(data_df)} rows x {ncol(data_df)} cols"
  )

  # Filter to valid classes only (exclude QC, NA, etc.)
  valid_classes <- acc_config$classes  # e.g. c(0, 1)
  if (!is.null(valid_classes)) {
    keep <- data_df$Class %in% valid_classes
    n_removed <- sum(!keep)
    data_df <- data_df[keep, , drop = FALSE]
    if (n_removed > 0) {
      cli::cli_alert_info(
        paste0(
          "Filtered to classes ",
          paste(valid_classes, collapse = ","),
          ": ", nrow(data_df), " kept, ",
          n_removed, " removed"
        )
      )
    }
  }

  # Extract feature columns (M1, M2, ...)
  feat_cols <- grep("^M[0-9]+$", names(data_df), value = TRUE)
  if (length(feat_cols) < 5) {
    cli::cli_alert_danger(
      "{acc_id}: too few feature columns ({length(feat_cols)})"
    )
    return(NULL)
  }

  X <- as.matrix(data_df[, feat_cols])
  X <- apply(X, 2, as.numeric)
  rownames(X) <- as.character(data_df$SampleID)

  # Get metabolite names from Peak sheet
  peak_df <- tryCatch(
    readxl::read_excel(xlsx_file, sheet = "Peak"),
    error = function(e) NULL
  )
  if (!is.null(peak_df) && "Name" %in% names(peak_df)) {
    feat_names <- make.unique(as.character(peak_df$Name))
    if (length(feat_names) == ncol(X)) {
      colnames(X) <- feat_names
    }
  }

  # Class labels — may be integer (0/1) or string
  class_vals <- data_df$Class
  if (is.numeric(class_vals) &&
      all(class_vals %in% c(0, 1), na.rm = TRUE)) {
    y_raw <- as.integer(class_vals)
    already_binary <- TRUE
  } else {
    # String classes — pass through to binarize_classes
    y_raw <- as.character(class_vals)
    already_binary <- FALSE
  }

  # Sample info: keep all non-feature columns
  meta_cols <- setdiff(names(data_df), feat_cols)
  sample_info <- as.data.frame(data_df[, meta_cols])

  feature_info <- data.frame(
    feature_id = colnames(X),
    original_label = feat_cols,
    stringsAsFactors = FALSE
  )

  cli::cli_alert_success(
    paste0(
      acc_id, ": ", nrow(X), " samples x ", ncol(X),
      " features, classes ",
      paste(names(table(y_raw)), collapse = "/")
    )
  )

  list(
    X = X,
    y = y_raw,
    sample_info = sample_info,
    feature_info = feature_info,
    y_already_binary = already_binary
  )
}

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

    # Look up config-specified class column for this accession
    acc_config <- NULL
    for (a in config$datasets$metabolights$accessions) {
      if (a$id == acc_id) { acc_config <- a; break }
    }

    # Determine which factor column to use
    if (!is.null(acc_config$class_column) && acc_config$class_column %in% names(samp)) {
      factor_cols <- acc_config$class_column
    } else {
      factor_cols <- get_matching_colnames(
        samp,
        "factor|group|class|disease|condition|diagnosis|status"
      )
    }

    if (length(factor_cols) > 0) {
      # Match sample names: MAF data column names are sample IDs
      sample_ids <- rownames(X)  # These are the MAF column names (= sample IDs)
      samp_name_col <- if ("Sample Name" %in% names(samp)) "Sample Name" else {
        # Try Source Name as fallback
        sn <- grep("^Source Name$|^Sample Name$|^Sample$", names(samp), value = TRUE)
        if (length(sn) > 0) sn[1] else NULL
      }

      if (!is.null(samp_name_col)) {
        matched <- match(sample_ids, as.character(samp[[samp_name_col]]))
        if (sum(!is.na(matched)) > 0) {
          sample_info <- samp[matched[!is.na(matched)], , drop = FALSE]
          y_raw <- samp[[factor_cols[1]]][matched]
          y <- as.character(y_raw)
          cli::cli_alert_info("Matched {sum(!is.na(matched))}/{length(sample_ids)} samples via '{samp_name_col}'")
        }
      }

      # Fallback: if no match and same number of rows, assume aligned
      if (is.null(y) && nrow(samp) == length(sample_ids)) {
        sample_info <- samp
        y <- as.character(samp[[factor_cols[1]]])
        cli::cli_alert_info("Assumed row-aligned sample sheet ({nrow(samp)} rows)")
      }

      # Apply sample filter if specified in config (e.g. keep only "biological material")
      if (!is.null(y) && !is.null(acc_config$filter_column) &&
          acc_config$filter_column %in% names(samp) && !is.null(samp_name_col)) {
        filter_vals <- samp[[acc_config$filter_column]]
        filter_match <- match(sample_ids, as.character(samp[[samp_name_col]]))
        keep <- rep(TRUE, length(sample_ids))
        keep[!is.na(filter_match)] <- filter_vals[filter_match[!is.na(filter_match)]] == acc_config$filter_value
        keep[is.na(keep)] <- FALSE
        if (sum(keep) < length(keep)) {
          X <- X[keep, , drop = FALSE]
          y <- y[keep]
          if (!is.null(sample_info)) {
            # re-align sample_info
            sample_info <- sample_info[seq_len(sum(keep)), , drop = FALSE]
          }
          cli::cli_alert_info("Filtered by '{acc_config$filter_column}' = '{acc_config$filter_value}': {sum(keep)}/{length(keep)} samples kept")
        }
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
#' MW REST API returns nested JSON:
#'   data: {row_id: {metabolite_name, DATA: {sample_id: value, ...}}, ...}
#'   factors: {row_id: {local_sample_id, factors: "Key:Val | Key:Val", ...}, ...}
#' @param acc_id Accession ID
#' @param raw_dir Raw data directory
#' @return List with X (matrix), y (vector), sample_info, feature_info
parse_mw <- function(acc_id, raw_dir) {
  study_dir <- file.path(raw_dir, acc_id)
  data_file <- file.path(study_dir, paste0(acc_id, "_data.json"))
  factors_file <- file.path(study_dir, paste0(acc_id, "_factors.json"))

  if (!file.exists(data_file)) {
    cli::cli_alert_danger("{acc_id}: no data file found")
    return(NULL)
  }

  # --- Parse data matrix (nested JSON) ----------------------------------------
  data_json <- jsonlite::fromJSON(data_file, simplifyDataFrame = FALSE)

  if (length(data_json) == 0) {
    cli::cli_alert_danger("{acc_id}: empty data JSON")
    return(NULL)
  }

  # Extract: each element has $metabolite_name and $DATA (named list sample->value)
  first_entry <- data_json[[1]]
  if (!is.null(first_entry$DATA)) {
    # Nested MW format: {id: {metabolite_name, DATA: {sample: val}}}
    sample_ids <- names(first_entry$DATA)
    metabolite_names <- character(length(data_json))
    mat <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(data_json))

    for (i in seq_along(data_json)) {
      entry <- data_json[[i]]
      metabolite_names[i] <- entry$metabolite_name %||% paste0("F", i)
      vals <- unlist(entry$DATA)
      # Ensure same sample order
      mat[, i] <- suppressWarnings(as.numeric(vals[sample_ids]))
    }

    # Trim leading/trailing whitespace from sample IDs (MW sometimes has spaces)
    sample_ids <- trimws(sample_ids)
    metabolite_names <- make.unique(metabolite_names)
    rownames(mat) <- sample_ids
    colnames(mat) <- metabolite_names
    X <- mat
    cli::cli_alert_info("Parsed MW nested JSON: {nrow(X)} samples × {ncol(X)} metabolites")
  } else if (is.data.frame(data_json)) {
    # Direct data frame format (rare)
    num_cols <- vapply(data_json, is.numeric, logical(1))
    if (sum(num_cols) < 5) {
      cli::cli_alert_danger("{acc_id}: too few numeric columns")
      return(NULL)
    }
    X <- as.matrix(data_json[, num_cols])
    sample_ids <- if ("sample_id" %in% names(data_json)) data_json$sample_id else paste0("S", seq_len(nrow(X)))
    rownames(X) <- sample_ids
  } else {
    cli::cli_alert_danger("{acc_id}: could not parse data JSON format")
    return(NULL)
  }

  if (ncol(X) < 5) {
    cli::cli_alert_danger("{acc_id}: too few features ({ncol(X)})")
    return(NULL)
  }

  # --- Parse factors (nested JSON with pipe-delimited factors string) ---------
  y <- NULL
  sample_info <- NULL

  if (file.exists(factors_file)) {
    factors_json <- tryCatch(
      jsonlite::fromJSON(factors_file, simplifyDataFrame = FALSE),
      error = function(e) NULL
    )

    if (!is.null(factors_json) && is.list(factors_json)) {
      # Build data.frame from nested list
      flist <- lapply(factors_json, function(entry) {
        sid <- trimws(entry$local_sample_id %||% "")
        fstring <- entry$factors %||% ""
        list(local_sample_id = sid, factors_raw = fstring)
      })
      fdf <- do.call(rbind, lapply(flist, as.data.frame, stringsAsFactors = FALSE))

      # Parse pipe-delimited factors into columns
      parsed_factors <- strsplit(fdf$factors_raw, " \\| ")
      all_keys <- unique(unlist(lapply(parsed_factors, function(x) sub(":.*", "", x))))

      for (k in all_keys) {
        fdf[[k]] <- vapply(parsed_factors, function(x) {
          hit <- grep(paste0("^", k, ":"), x, value = TRUE)
          if (length(hit) > 0) sub(paste0("^", k, ":"), "", hit[1]) else NA_character_
        }, character(1))
      }

      # Match to X sample IDs
      matched <- match(rownames(X), fdf$local_sample_id)
      n_matched <- sum(!is.na(matched))
      cli::cli_alert_info("Factor matching: {n_matched}/{nrow(X)} samples")

      if (n_matched > 0) {
        sample_info <- fdf[matched[!is.na(matched)], , drop = FALSE]

        # Determine class column from config or auto-detect
        acc_config <- NULL
        for (a in config$datasets$metabolomics_workbench$accessions) {
          if (a$id == acc_id) { acc_config <- a; break }
        }

        class_col <- NULL
        if (!is.null(acc_config$class_factor) && acc_config$class_factor %in% names(fdf)) {
          class_col <- acc_config$class_factor
        } else {
          # Auto-detect: look for columns with 2-4 unique values
          for (k in all_keys) {
            uv <- unique(fdf[[k]][!is.na(fdf[[k]])])
            if (length(uv) >= 2 && length(uv) <= 4) { class_col <- k; break }
          }
        }

        if (!is.null(class_col)) {
          y_full <- rep(NA_character_, nrow(X))
          y_full[!is.na(matched)] <- fdf[[class_col]][matched[!is.na(matched)]]
          y <- y_full
          cli::cli_alert_info("Class column: '{class_col}' → {paste(names(table(y)), collapse=', ')}")
        }
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
      # Log2(x + 1) transform + autoscaling (mean-center, unit variance)
      X_log <- log2(X + 1)
      scale(X_log, center = TRUE, scale = TRUE)
    },

    "log_pareto" = {
      # Log2(x + 1) transform + Pareto scaling (divide by sqrt(sd))
      X_log <- log2(X + 1)
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
      X_log <- log2(X_pqn + 1)
      scale(X_log, center = TRUE, scale = TRUE)
    },

    # Default fallback
    {
      cli::cli_alert_warning("Unknown method '{method}', using log_auto")
      X_log <- log2(X + 1)
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
    if (info$source == "CIMCB") {
      # Find config entry for this accession
      acc_cfg <- NULL
      for (a in config$datasets$cimcb$accessions) {
        if (a$id == acc_id) { acc_cfg <- a; break }
      }
      parse_cimcb(acc_id, raw_dir, acc_cfg)
    } else if (info$source == "MetaboLights") {
      parse_metabolights(acc_id, raw_dir)
    } else {
      parse_mw(acc_id, raw_dir)
    }
  }, error = function(e) {
    cli::cli_alert_danger("Parse error: {e$message}")
    NULL
  })

  if (is.null(parsed)) {
    cli::cli_alert_warning(
      "{acc_id}: parsing failed, skipping"
    )
    datasets_summary[[acc_id]] <- list(
      id = acc_id, status = "parse_failed"
    )
    next
  }

  # Apply QC filters
  filtered <- apply_qc_filters(parsed, config)
  if (is.null(filtered)) {
    cli::cli_alert_warning("{acc_id}: failed QC, skipping")
    datasets_summary[[acc_id]] <- list(
      id = acc_id, status = "qc_failed"
    )
    next
  }

  # Binarize classes (skip for CIMCB — already binary)
  if (isTRUE(filtered$y_already_binary)) {
    y <- as.integer(filtered$y)
    valid_idx <- !is.na(y)
    X <- filtered$X[valid_idx, , drop = FALSE]
    y <- y[valid_idx]
    class_map <- c(control = "0", case = "1")
    sample_info <- filtered$sample_info
    if (!is.null(sample_info) &&
        nrow(sample_info) == length(valid_idx)) {
      sample_info <- sample_info[valid_idx, , drop = FALSE]
    }
  } else {
    class_result <- binarize_classes(filtered$y)
    if (!class_result$valid) {
      cli::cli_alert_warning(
        "{acc_id}: cannot binarize classes, skipping"
      )
      datasets_summary[[acc_id]] <- list(
        id = acc_id, status = "class_failed"
      )
      next
    }
    valid_idx <- !is.na(class_result$y_binary)
    X <- filtered$X[valid_idx, , drop = FALSE]
    y <- class_result$y_binary[valid_idx]
    class_map <- class_result$class_map
    sample_info <- filtered$sample_info
    if (!is.null(sample_info) &&
        nrow(sample_info) == length(valid_idx)) {
      sample_info <- sample_info[valid_idx, , drop = FALSE]
    }
  }

  # Check minimum samples per group
  tab <- table(y)
  min_spg <- config$acceptance_criteria$min_samples_per_group
  if (any(tab < min_spg)) {
    cli::cli_alert_warning(
      paste0(
        acc_id, ": insufficient samples/group (",
        paste(tab, collapse = "/"), "), skipping"
      )
    )
    datasets_summary[[acc_id]] <- list(
      id = acc_id, status = "too_few_samples"
    )
    next
  }

  # Impute
  X_pre_impute <- X
  X <- impute_missing(X)

  # Normalize
  X_norm <- normalize_matrix(
    X, method = config$simulation$base$preprocessing
  )

  # Final validation
  validate_data(X_norm, y, label = acc_id)

  # Save
  processed <- list(
    X = X_norm,
    X_raw = X_pre_impute,
    X_imputed = X,
    y = y,
    class_map = class_map,
    feature_info = filtered$feature_info,
    sample_info = sample_info,
    source = info$source,
    platform = info$platform,
    accession = acc_id,
    description = info$description
  )

  out_file <- file.path(
    proc_dir, paste0(acc_id, "_processed.rds")
  )
  save_result(
    processed, out_file,
    metadata = list(
      n = nrow(X_norm), p = ncol(X_norm),
      n_case = sum(y == 1), n_control = sum(y == 0)
    )
  )

  datasets_summary[[acc_id]] <- list(
    id = acc_id, status = "success",
    source = info$source, platform = info$platform,
    n = nrow(X_norm), p = ncol(X_norm),
    n_case = sum(y == 1), n_control = sum(y == 0)
  )

  cli::cli_alert_success(
    "{acc_id}: {nrow(X_norm)} samples x {ncol(X_norm)} features"
  )
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
