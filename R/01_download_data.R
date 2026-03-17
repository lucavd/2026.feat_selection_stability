# ==============================================================================
# Script 01 — Download public metabolomics datasets
# ==============================================================================
# Purpose: Download feature tables from MetaboLights and Metabolomics Workbench.
# Input:   config/config.yaml (dataset accessions)
# Output:  data/raw/<accession>/ directories with raw data files
#          data/raw/download_manifest.rds (status log)
# ==============================================================================

cat("=== Script 01: Download Data ===\n")

suppressPackageStartupMessages({
  library(here)
  library(cli)
  library(httr)
  library(jsonlite)
  library(curl)
})

source(here("R", "utils", "helpers.R"))
config <- load_config()
setup_logging("Script 01: Download Data", config)

raw_dir <- get_path(config, "raw_data")

# --- Download manifest --------------------------------------------------------
# Tracks what was downloaded, when, and status
manifest <- list()

# ==============================================================================
# MetaboLights downloads
# ==============================================================================

cli::cli_h1("MetaboLights datasets")

download_metabolights <- function(accession_info, config) {
  acc_id   <- accession_info$id
  platform <- accession_info$platform
  desc     <- accession_info$description

  cli::cli_h2("{acc_id} — {desc} ({platform})")

  out_dir <- file.path(raw_dir, acc_id)
  fs::dir_create(out_dir)

  result <- list(
    id = acc_id, source = "MetaboLights", platform = platform,
    description = desc, status = "pending",
    timestamp = Sys.time(), files = character(0)
  )

  tryCatch({
    # Step 1: Get study metadata via REST API
    api_base <- config$datasets$metabolights$api_base
    study_url <- paste0(api_base, "/studies/", acc_id)

    cli::cli_alert_info("Fetching study metadata from API...")
    resp <- httr::GET(
      study_url,
      httr::timeout(config$datasets$download$timeout_seconds),
      httr::add_headers(Accept = "application/json")
    )

    if (httr::status_code(resp) != 200) {
      cli::cli_alert_warning("API returned status {httr::status_code(resp)}, trying FTP...")
    } else {
      meta <- httr::content(resp, as = "text", encoding = "UTF-8")
      writeLines(meta, file.path(out_dir, "study_metadata.json"))
      cli::cli_alert_success("Study metadata saved")
    }

    # Step 2: List study files
    files_url <- paste0(api_base, "/studies/", acc_id, "/files?include_raw_data=false")
    resp_files <- httr::GET(
      files_url,
      httr::timeout(config$datasets$download$timeout_seconds)
    )

    target_files <- c()
    if (httr::status_code(resp_files) == 200) {
      files_info <- httr::content(resp_files, as = "parsed")
      # Look for metabolite assignment files (MAF) and sample sheets
      if (!is.null(files_info$study)) {
        all_files <- unlist(lapply(files_info$study, function(x) x$file))
      } else {
        all_files <- character(0)
      }

      # Priority: MAF files (m_*.tsv), sample sheet (s_*.txt), assay (a_*.txt)
      maf_files <- grep("^m_.*\\.tsv$", all_files, value = TRUE)
      sample_files <- grep("^s_.*\\.txt$", all_files, value = TRUE)
      assay_files <- grep("^a_.*\\.txt$", all_files, value = TRUE)
      target_files <- c(maf_files, sample_files, assay_files)
      cli::cli_alert_info("Found {length(target_files)} target files")
    }

    # Step 3: Download via FTP
    ftp_base <- config$datasets$metabolights$ftp_base
    ftp_study <- paste0(ftp_base, "/", acc_id)

    # If no specific files found, try known patterns
    if (length(target_files) == 0) {
      cli::cli_alert_info("No file list from API; trying known file patterns via FTP")
      target_files <- c(
        paste0("m_", acc_id, "_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv"),
        paste0("m_", acc_id, "_metabolite_profiling_mass_spectrometry_v2_maf.tsv"),
        paste0("s_", acc_id, ".txt"),
        paste0("a_", acc_id, ".txt")
      )
    }

    downloaded <- character(0)
    for (fname in target_files) {
      ftp_url <- paste0(ftp_study, "/", fname)
      local_file <- file.path(out_dir, fname)

      for (attempt in seq_len(config$datasets$download$max_retries)) {
        dl_ok <- tryCatch({
          curl::curl_download(ftp_url, local_file,
                              quiet = TRUE,
                              handle = curl::new_handle(
                                connecttimeout = config$datasets$download$timeout_seconds
                              ))
          TRUE
        }, error = function(e) {
          if (attempt < config$datasets$download$max_retries) {
            Sys.sleep(config$datasets$download$retry_backoff_seconds * attempt)
          }
          FALSE
        })
        if (dl_ok) break
      }

      if (dl_ok && file.exists(local_file) && file.size(local_file) > 0) {
        downloaded <- c(downloaded, fname)
        cli::cli_alert_success("  Downloaded: {fname}")
      }
    }

    # Step 4: Try alternative ISA-Tab download if nothing found
    if (length(downloaded) == 0) {
      cli::cli_alert_info("Trying ISA-Tab bundle download...")
      isa_url <- paste0(api_base, "/studies/", acc_id, "/download")
      isa_file <- file.path(out_dir, paste0(acc_id, "_ISATab.zip"))
      dl_ok <- tryCatch({
        httr::GET(isa_url,
                  httr::write_disk(isa_file, overwrite = TRUE),
                  httr::timeout(config$datasets$download$timeout_seconds))
        if (file.exists(isa_file) && file.size(isa_file) > 100) {
          utils::unzip(isa_file, exdir = out_dir)
          cli::cli_alert_success("ISA-Tab bundle extracted")
          downloaded <- list.files(out_dir, pattern = "\\.(tsv|txt)$")
          TRUE
        } else {
          FALSE
        }
      }, error = function(e) FALSE)
    }

    result$files  <- downloaded
    result$status <- if (length(downloaded) > 0) "success" else "no_files"
    if (result$status == "success") {
      cli::cli_alert_success("{acc_id}: {length(downloaded)} files downloaded")
    } else {
      cli::cli_alert_warning("{acc_id}: no data files downloaded")
    }

  }, error = function(e) {
    result$status  <<- "error"
    result$message <<- e$message
    cli::cli_alert_danger("{acc_id}: {e$message}")
  })

  result
}

# Download each MetaboLights dataset
for (acc_info in config$datasets$metabolights$accessions) {
  manifest[[acc_info$id]] <- download_metabolights(acc_info, config)
}

# ==============================================================================
# Metabolomics Workbench downloads
# ==============================================================================

cli::cli_h1("Metabolomics Workbench datasets")

download_mw <- function(accession_info, config) {
  acc_id   <- accession_info$id
  platform <- accession_info$platform
  desc     <- accession_info$description

  cli::cli_h2("{acc_id} — {desc} ({platform})")

  out_dir <- file.path(raw_dir, acc_id)
  fs::dir_create(out_dir)

  result <- list(
    id = acc_id, source = "MetabolomicsWorkbench", platform = platform,
    description = desc, status = "pending",
    timestamp = Sys.time(), files = character(0)
  )

  tryCatch({
    api_base <- config$datasets$metabolomics_workbench$api_base

    # Step 1: Get study summary
    summary_url <- paste0(api_base, "/study/study_id/", acc_id, "/summary")
    cli::cli_alert_info("Fetching study summary...")
    resp <- httr::GET(summary_url,
                      httr::timeout(config$datasets$download$timeout_seconds))

    if (httr::status_code(resp) == 200) {
      summary_text <- httr::content(resp, as = "text", encoding = "UTF-8")
      writeLines(summary_text, file.path(out_dir, "study_summary.json"))
    }

    # Step 2: Download data matrix (metabolite concentrations)
    # MW REST API: study_id/ST.../datatable/...
    data_url <- paste0(api_base, "/study/study_id/", acc_id,
                       "/data")
    cli::cli_alert_info("Fetching data matrix...")
    resp_data <- httr::GET(data_url,
                           httr::timeout(config$datasets$download$timeout_seconds))

    downloaded <- character(0)
    if (httr::status_code(resp_data) == 200) {
      data_text <- httr::content(resp_data, as = "text", encoding = "UTF-8")
      if (nchar(data_text) > 100) {
        data_file <- file.path(out_dir, paste0(acc_id, "_data.json"))
        writeLines(data_text, data_file)
        downloaded <- c(downloaded, basename(data_file))
        cli::cli_alert_success("Data matrix saved")
      }
    }

    # Step 3: Get sample factors (phenotype/class info)
    factors_url <- paste0(api_base, "/study/study_id/", acc_id, "/factors")
    resp_factors <- httr::GET(factors_url,
                              httr::timeout(config$datasets$download$timeout_seconds))

    if (httr::status_code(resp_factors) == 200) {
      factors_text <- httr::content(resp_factors, as = "text", encoding = "UTF-8")
      if (nchar(factors_text) > 10) {
        factors_file <- file.path(out_dir, paste0(acc_id, "_factors.json"))
        writeLines(factors_text, factors_file)
        downloaded <- c(downloaded, basename(factors_file))
        cli::cli_alert_success("Sample factors saved")
      }
    }

    # Step 4: Try metabolite data table
    metab_url <- paste0(api_base, "/study/study_id/", acc_id,
                        "/untarg_metabolite_data")
    resp_metab <- httr::GET(metab_url,
                            httr::timeout(config$datasets$download$timeout_seconds))

    if (httr::status_code(resp_metab) == 200) {
      metab_text <- httr::content(resp_metab, as = "text", encoding = "UTF-8")
      if (nchar(metab_text) > 100) {
        metab_file <- file.path(out_dir, paste0(acc_id, "_metabolites.json"))
        writeLines(metab_text, metab_file)
        downloaded <- c(downloaded, basename(metab_file))
        cli::cli_alert_success("Metabolite data saved")
      }
    }

    result$files  <- downloaded
    result$status <- if (length(downloaded) > 0) "success" else "no_files"
    if (result$status == "success") {
      cli::cli_alert_success("{acc_id}: {length(downloaded)} files downloaded")
    } else {
      cli::cli_alert_warning("{acc_id}: no data files downloaded")
    }

  }, error = function(e) {
    result$status  <<- "error"
    result$message <<- e$message
    cli::cli_alert_danger("{acc_id}: {e$message}")
  })

  result
}

# Download each MW dataset
for (acc_info in config$datasets$metabolomics_workbench$accessions) {
  manifest[[acc_info$id]] <- download_mw(acc_info, config)
}

# ==============================================================================
# Save manifest
# ==============================================================================

cli::cli_h1("Download Summary")

manifest_df <- do.call(rbind, lapply(manifest, function(m) {
  data.frame(
    id = m$id, source = m$source, platform = m$platform,
    status = m$status, n_files = length(m$files),
    stringsAsFactors = FALSE
  )
}))

print(manifest_df)

save_result(
  data = manifest,
  path = file.path(raw_dir, "download_manifest.rds"),
  metadata = list(script = "01_download_data.R",
                  n_datasets = length(manifest),
                  n_success = sum(manifest_df$status == "success"))
)

# --- Session info -------------------------------------------------------------
finalize_logging("Script 01", session_file = "results/session_info_01.txt")

cat("\n=== Script 01: Complete ===\n")
