# ==============================================================================
# Utility functions — stability_metrics.R
# ==============================================================================
# Wrappers for computing stability metrics on feature selection results.
# Primary metric: Nogueira index (Nogueira et al. 2018, JMLR 18(174):1-54)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(cli)
})

# --- Nogueira stability index ------------------------------------------------

#' Compute Nogueira stability index with confidence interval
#'
#' Implements the corrected stability measure from:
#' Nogueira, Sechidis & Brown (2018). On the stability of feature selection
#' algorithms. JMLR 18(174):1-54.
#'
#' @param selection_matrix Binary matrix (B × p): rows = bootstrap runs,
#'   cols = features, 1 = selected, 0 = not selected.
#' @param alpha Significance level for CI (default 0.05 → 95% CI)
#' @return List with: estimate, variance, ci_lower, ci_upper, n_bootstrap, p
compute_nogueira <- function(selection_matrix, alpha = 0.05) {
  stopifnot(is.matrix(selection_matrix))
  stopifnot(all(selection_matrix %in% c(0, 1)))

  B <- nrow(selection_matrix)  # Number of bootstrap samples

p <- ncol(selection_matrix)  # Number of features

  if (B < 2) {
    cli::cli_alert_warning("Nogueira index requires B >= 2; returning NA")
    return(list(estimate = NA_real_, variance = NA_real_,
                ci_lower = NA_real_, ci_upper = NA_real_,
                n_bootstrap = B, p = p))
  }

  # Feature selection frequencies
  pf <- colMeans(selection_matrix)          # p_hat_j for each feature
  kbar <- mean(rowSums(selection_matrix))   # Average number selected

  # Nogueira index (Eq. 4 in paper)
  # phi = 1 - (1 / (kbar/p * (1 - kbar/p))) * (1/p) * sum(pf * (1 - pf)) * (B / (B-1))
  kbar_ratio <- kbar / p

  # Handle degenerate cases
  denom <- kbar_ratio * (1 - kbar_ratio)
  if (denom < 1e-10) {
    # All features always selected or never selected → perfect stability
    # but trivially so
    return(list(estimate = 1.0, variance = 0.0,
                ci_lower = 1.0, ci_upper = 1.0,
                n_bootstrap = B, p = p))
  }

  phi_hat <- 1 - (1 / denom) * (1 / p) * sum(pf * (1 - pf)) * (B / (B - 1))

  if (B == 2) {
    return(list(estimate = phi_hat, variance = NA_real_,
                ci_lower = NA_real_, ci_upper = NA_real_,
                n_bootstrap = B, p = p))
  }

  # Variance of the estimator (Theorem 7 in paper)
  # Based on the leave-one-out jackknife variance
  phi_jack <- numeric(B)
  for (b in seq_len(B)) {
    sm_loo <- selection_matrix[-b, , drop = FALSE]
    pf_loo <- colMeans(sm_loo)
    kbar_loo <- mean(rowSums(sm_loo))
    kbar_ratio_loo <- kbar_loo / p
    denom_loo <- kbar_ratio_loo * (1 - kbar_ratio_loo)
    if (denom_loo < 1e-10) {
      phi_jack[b] <- 1.0
    } else {
      phi_jack[b] <- 1 - (1 / denom_loo) * (1 / p) *
        sum(pf_loo * (1 - pf_loo)) * ((B - 1) / (B - 2))
    }
  }

  var_hat <- ((B - 1) / B) * sum((phi_jack - mean(phi_jack))^2)
  z <- qnorm(1 - alpha / 2)

  list(
    estimate   = phi_hat,
    variance   = var_hat,
    ci_lower   = phi_hat - z * sqrt(var_hat),
    ci_upper   = phi_hat + z * sqrt(var_hat),
    n_bootstrap = B,
    p          = p
  )
}

# --- Other stability metrics (via stabm where available) ----------------------

#' Compute Jaccard stability index
#'
#' @param selection_matrix Binary matrix (B × p)
#' @return Numeric mean pairwise Jaccard similarity
compute_jaccard <- function(selection_matrix) {
  B <- nrow(selection_matrix)
  if (B < 2) return(NA_real_)

  # Convert rows to sets (indices of selected features)
  sets <- lapply(seq_len(B), function(b) which(selection_matrix[b, ] == 1))

  # Pairwise Jaccard
  pairs <- combn(B, 2)
  jaccards <- apply(pairs, 2, function(idx) {
    a <- sets[[idx[1]]]
    b <- sets[[idx[2]]]
    inter <- length(intersect(a, b))
    union <- length(union(a, b))
    if (union == 0) return(1.0)  # Both empty
    inter / union
  })

  mean(jaccards)
}

#' Compute Kuncheva index
#'
#' Kuncheva (2007) consistency index, corrected for chance.
#'
#' @param selection_matrix Binary matrix (B × p)
#' @return Numeric Kuncheva index
compute_kuncheva <- function(selection_matrix) {
  B <- nrow(selection_matrix)
  p <- ncol(selection_matrix)
  if (B < 2) return(NA_real_)

  k_sizes <- rowSums(selection_matrix)

  # Pairwise
  pairs <- combn(B, 2)
  kuncheva_vals <- apply(pairs, 2, function(idx) {
    s1 <- which(selection_matrix[idx[1], ] == 1)
    s2 <- which(selection_matrix[idx[2], ] == 1)
    r <- length(intersect(s1, s2))
    k1 <- length(s1)
    k2 <- length(s2)
    denom <- k1 * k2 / p
    if (abs(min(k1, k2) - denom) < 1e-10) return(NA_real_)
    (r - denom) / (min(k1, k2) - denom)
  })

  mean(kuncheva_vals, na.rm = TRUE)
}

#' Compute Dice similarity
#'
#' @param selection_matrix Binary matrix (B × p)
#' @return Numeric mean pairwise Dice coefficient
compute_dice <- function(selection_matrix) {
  B <- nrow(selection_matrix)
  if (B < 2) return(NA_real_)

  sets <- lapply(seq_len(B), function(b) which(selection_matrix[b, ] == 1))

  pairs <- combn(B, 2)
  dices <- apply(pairs, 2, function(idx) {
    a <- sets[[idx[1]]]
    b <- sets[[idx[2]]]
    inter <- length(intersect(a, b))
    total <- length(a) + length(b)
    if (total == 0) return(1.0)
    2 * inter / total
  })

  mean(dices)
}

#' Compute Spearman rank correlation of feature rankings
#'
#' For methods that produce importance scores, not just binary selections.
#'
#' @param importance_matrix Numeric matrix (B × p) of importance scores
#' @return Numeric mean pairwise Spearman correlation
compute_spearman_rank <- function(importance_matrix) {
  B <- nrow(importance_matrix)
  if (B < 2) return(NA_real_)

  # Rank each row
  rank_matrix <- t(apply(importance_matrix, 1, rank, ties.method = "average"))

  pairs <- combn(B, 2)
  cors <- apply(pairs, 2, function(idx) {
    cor(rank_matrix[idx[1], ], rank_matrix[idx[2], ], method = "spearman")
  })

  mean(cors, na.rm = TRUE)
}

# --- Master function ----------------------------------------------------------

#' Compute all stability metrics
#'
#' @param selection_matrix Binary matrix (B × p)
#' @param importance_matrix Optional numeric matrix (B × p) for rank-based metrics
#' @param metrics_config Config list for which metrics to compute
#' @return Named list of metric results
compute_all_stability <- function(selection_matrix,
                                  importance_matrix = NULL,
                                  metrics_config = NULL) {
  results <- list()

  # Primary metrics
  results$nogueira <- compute_nogueira(selection_matrix)
  results$jaccard  <- compute_jaccard(selection_matrix)

  # Secondary metrics
  results$kuncheva <- compute_kuncheva(selection_matrix)
  results$dice     <- compute_dice(selection_matrix)

  if (!is.null(importance_matrix)) {
    results$spearman_rank <- compute_spearman_rank(importance_matrix)
  }

  # Summary stats
  results$n_bootstrap <- nrow(selection_matrix)
  results$p           <- ncol(selection_matrix)
  results$mean_n_selected <- mean(rowSums(selection_matrix))
  results$sd_n_selected   <- sd(rowSums(selection_matrix))

  results
}
