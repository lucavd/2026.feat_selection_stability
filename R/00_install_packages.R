# ==============================================================================
# Script 00 — Install and verify all required packages
# ==============================================================================
# Purpose: Ensure all dependencies are installed before running the pipeline.
# Usage:   Rscript R/00_install_packages.R
# ==============================================================================

cat("=== Script 00: Installing packages ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- CRAN packages -----------------------------------------------------------
cran_packages <- c(

  # Core infrastructure
  "here", "yaml", "cli", "fs", "withr",
  # Data manipulation

  "data.table", "dplyr", "tidyr", "purrr", "tibble", "readr", "stringr",
  # Statistical / ML

  "glmnet", "Boruta", "ranger", "stabs", "knockoff",
  "horseshoe", "spikeslab",
  "xgboost",
  # Simulation

  "mvtnorm", "corpcor", "MASS",
  # Parallelization

  "future", "furrr", "progressr",
  # Metrics
  "stabm", "pROC", "caret",
  # Visualization
  "ggplot2", "patchwork", "scales", "viridis", "ggrepel",
  "ComplexHeatmap",
  # SHAP
  "shapviz", "treeshap",
  # Imputation
  "missForest", "VIM",
  # Reporting
  "knitr", "rmarkdown",
  # API / download
  "httr", "jsonlite", "curl"
)

# --- Bioconductor packages ---------------------------------------------------
bioc_packages <- c(
  "metabolomicsWorkbenchR",
  "SummarizedExperiment",
  "Biobase"
)

# --- GitHub packages ----------------------------------------------------------
github_packages <- list(
  # metabolighteR for MetaboLights API
  list(repo = "aberHRML/metabolighteR", pkg = "metabolighteR")
)

# --- Install CRAN packages ----------------------------------------------------
install_cran <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_alert_info("Installing CRAN package: {pkg}")
      tryCatch(
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         quiet = TRUE, dependencies = TRUE),
        error = function(e) {
          cli::cli_alert_danger("Failed to install {pkg}: {e$message}")
        }
      )
    }
  }
}

# --- Install Bioconductor packages -------------------------------------------
install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_alert_info("Installing Bioconductor package: {pkg}")
      tryCatch(
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE),
        error = function(e) {
          cli::cli_alert_danger("Failed to install {pkg}: {e$message}")
        }
      )
    }
  }
}

# --- Install GitHub packages --------------------------------------------------
install_gh <- function(pkgs) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
  for (entry in pkgs) {
    if (!requireNamespace(entry$pkg, quietly = TRUE)) {
      cli::cli_alert_info("Installing GitHub package: {entry$repo}")
      tryCatch(
        remotes::install_github(entry$repo, quiet = TRUE, upgrade = "never"),
        error = function(e) {
          cli::cli_alert_danger("Failed to install {entry$repo}: {e$message}")
        }
      )
    }
  }
}

# --- Run installations --------------------------------------------------------
cli::cli_h1("Installing CRAN packages")
install_cran(cran_packages)

cli::cli_h1("Installing Bioconductor packages")
install_bioc(bioc_packages)

cli::cli_h1("Installing GitHub packages")
install_gh(github_packages)

# --- Verify all packages ------------------------------------------------------
cli::cli_h1("Verification")

all_pkgs <- c(cran_packages, bioc_packages,
              vapply(github_packages, \(x) x$pkg, character(1)))

status <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
n_ok   <- sum(status)
n_fail <- sum(!status)

cli::cli_alert_success("{n_ok}/{length(all_pkgs)} packages available")

if (n_fail > 0) {
  cli::cli_alert_danger("Missing packages: {paste(names(status[!status]), collapse = ', ')}")
  cli::cli_alert_warning("Some methods may not run. Check installation logs above.")
} else {
  cli::cli_alert_success("All packages installed successfully.")
}

# --- Session info -------------------------------------------------------------
cli::cli_h1("Session Info")
writeLines(capture.output(sessionInfo()),
           here::here("results", "session_info_00.txt"))
cli::cli_alert_success("Session info saved to results/session_info_00.txt")

cat("\n=== Script 00: Complete ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
