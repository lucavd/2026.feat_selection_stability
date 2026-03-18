# ==============================================================================
# Utility functions — fs_methods.R
# ==============================================================================
# Wrapper functions for all 12 feature selection methods.
# Each wrapper returns a named list:
#   - selected:   logical vector (length p), TRUE if feature selected
#   - importance: numeric vector (length p), importance scores (NA if N/A)
#   - time:       elapsed time in seconds
#   - converged:  logical, TRUE if method ran successfully
#   - message:    character, error/warning message if any
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(cli)
})

# ==============================================================================
# 1. FILTER METHODS
# ==============================================================================

#' Wilcoxon rank-sum test with FDR correction
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with fdr_threshold
#' @return Standard FS result list
fs_wilcoxon_fdr <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    p_vals <- apply(X, 2, function(col) {
      tryCatch(
        wilcox.test(col[y == 1], col[y == 0])$p.value,
        error = function(e) NA_real_
      )
    })
    p_adj <- p.adjust(p_vals, method = "BH")
    selected <- !is.na(p_adj) & p_adj < params$fdr_threshold
    importance <- -log10(pmax(p_adj, 1e-300))  # Higher = more significant

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

#' Fold change filter
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with fc_threshold
#' @return Standard FS result list
fs_fold_change <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    mean_case    <- colMeans(X[y == 1, , drop = FALSE], na.rm = TRUE)
    mean_control <- colMeans(X[y == 0, , drop = FALSE], na.rm = TRUE)
    fc <- mean_case / pmax(mean_control, 1e-10)
    abs_log2fc <- abs(log2(fc))
    selected <- abs_log2fc >= log2(params$fc_threshold)
    importance <- abs_log2fc

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

#' Volcano filter (Wilcoxon FDR + fold change)
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with fdr_threshold, fc_threshold
#' @return Standard FS result list
fs_volcano <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    wilcox_res <- fs_wilcoxon_fdr(X, y, params)
    fc_res     <- fs_fold_change(X, y, params)
    selected <- wilcox_res$selected & fc_res$selected
    importance <- wilcox_res$importance * fc_res$importance

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

# ==============================================================================
# 2. EMBEDDED METHODS
# ==============================================================================

#' LASSO (L1-penalized logistic regression)
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with alpha, lambda_rule, nfolds
#' @return Standard FS result list
fs_lasso <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    fit <- glmnet::cv.glmnet(
      x = X, y = y, family = "binomial",
      alpha = params$alpha,
      nfolds = params$nfolds,
      type.measure = "deviance"
    )
    coefs <- as.numeric(coef(fit, s = params$lambda_rule))[-1]  # Drop intercept
    selected <- coefs != 0
    importance <- abs(coefs)

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

#' Elastic Net (alpha = 0.5 by default)
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with alpha, lambda_rule, nfolds
#' @return Standard FS result list
fs_elastic_net <- function(X, y, params) {
  # Same implementation as LASSO, different alpha
  fs_lasso(X, y, params)
}

#' Knockoff filter
#'
#' @param X Numeric matrix (n × p), requires n > p
#' @param y Binary response (0/1)
#' @param params List with fdr, statistic
#' @return Standard FS result list
fs_knockoff <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    # Knockoff requires n > p; check
    if (nrow(X) <= ncol(X)) {
      cli::cli_alert_warning("Knockoff: n <= p ({nrow(X)} <= {ncol(X)}), skipping")
      return(list(
        selected = rep(FALSE, ncol(X)),
        importance = rep(NA_real_, ncol(X)),
        time = as.numeric(proc.time()["elapsed"] - t0),
        converged = FALSE,
        message = "n <= p: knockoff requires n > p"
      ))
    }

    # Map statistic name to function
    stat_fn <- switch(params$statistic,
      "stat.glmnet_coefdiff" = knockoff::stat.glmnet_coefdiff,
      "stat.glmnet_lambdadiff" = knockoff::stat.glmnet_lambdadiff,
      knockoff::stat.glmnet_coefdiff  # default
    )

    result <- knockoff::knockoff.filter(
      X = X, y = y,
      fdr = params$fdr,
      statistic = stat_fn
    )

    selected <- seq_len(ncol(X)) %in% result$selected
    importance <- rep(0, ncol(X))
    importance[result$selected] <- 1

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

# ==============================================================================
# 3. WRAPPER METHODS
# ==============================================================================

#' Boruta feature selection
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (factor)
#' @param params List with maxRuns, pValue
#' @return Standard FS result list
fs_boruta <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    y_factor <- factor(y, levels = c(0, 1), labels = c("control", "case"))
    fit <- Boruta::Boruta(
      x = X, y = y_factor,
      maxRuns = params$maxRuns,
      pValue = params$pValue,
      doTrace = 0,
      num.threads = 1
    )
    # Resolve tentative features
    fit_final <- Boruta::TentativeRoughFix(fit)
    decision <- fit_final$finalDecision
    selected <- decision == "Confirmed"
    # Use median importance (ImpHistory)
    imp_history <- fit$ImpHistory
    importance <- apply(imp_history[, seq_len(ncol(X)), drop = FALSE], 2, median,
                        na.rm = TRUE)

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

#' Random Forest variable importance (ranger)
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (factor)
#' @param params List with num_trees, importance, top_k_method, top_n
#' @return Standard FS result list
fs_rf_importance <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    y_factor <- factor(y, levels = c(0, 1), labels = c("control", "case"))
    df <- data.frame(y = y_factor, X)

    fit <- ranger::ranger(
      y ~ ., data = df,
      num.trees = params$num_trees,
      importance = params$importance,
      probability = TRUE,
      num.threads = 1
    )

    importance <- fit$variable.importance

    # Select top features by elbow or fixed top_n
    selected <- select_top_k(importance, params$top_k_method, params$top_n)

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

# ==============================================================================
# 4. META METHODS
# ==============================================================================

#' Stability selection (via stabs package)
#'
#' Uses LASSO as the base learner.
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with cutoff, PFER, sampling_type, B
#' @return Standard FS result list
fs_stability_selection <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    fit <- stabs::stabsel(
      x = X, y = y,
      fitfun = stabs::glmnet.lasso,
      cutoff = params$cutoff,
      PFER = params$PFER,
      sampling.type = params$sampling_type,
      B = params$B
    )

    selected <- rep(FALSE, ncol(X))
    names(selected) <- colnames(X)
    selected[fit$selected] <- TRUE
    importance <- fit$max  # Maximum selection probability per feature

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

# ==============================================================================
# 5. BAYESIAN METHODS
# ==============================================================================

#' Horseshoe prior (continuous shrinkage)
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with method_tau, method_sigma, burn, nmc
#' @return Standard FS result list
fs_horseshoe <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    fit <- horseshoe::horseshoe(
      y = y, X = X,
      method.tau = params$method_tau,
      method.sigma = params$method_sigma,
      burn = params$burn,
      nmc = params$nmc
    )

    # Use posterior median of coefficients
    beta_hat <- fit$BetaHat
    # Selection: credible interval excludes zero
    ci_lower <- apply(fit$BetaSamples, 1, quantile, probs = 0.025)
    ci_upper <- apply(fit$BetaSamples, 1, quantile, probs = 0.975)
    selected <- (ci_lower > 0 | ci_upper < 0)
    importance <- abs(beta_hat)

    list(selected = selected, importance = as.numeric(importance),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

#' Spike-and-slab prior
#'
#' @param X Numeric matrix (n × p)
#' @param y Numeric response
#' @param params List with bigp_smalln, max_var
#' @return Standard FS result list
fs_spike_slab <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    fit <- spikeslab::spikeslab(
      formula = y ~ .,
      data = data.frame(y = y, X),
      bigp.smalln = params$bigp_smalln,
      max.var = min(params$max_var, ncol(X))
    )

    # gnet coefficients — non-zero means selected
    beta <- fit$gnet
    selected <- beta != 0
    importance <- abs(beta)

    list(selected = selected, importance = as.numeric(importance),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = "")
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

# ==============================================================================
# 6. ML METHODS
# ==============================================================================

#' SHAP-based selection via XGBoost
#'
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params List with max_depth, nrounds, eta, top_k_method, top_n
#' @return Standard FS result list
fs_shap_xgboost <- function(X, y, params) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

    # Use GPU if available, fall back to CPU (cached after first check)
    if (!exists(".xgb_device", envir = .GlobalEnv)) {
      device <- tryCatch({
        test_fit <- xgboost::xgb.train(
          list(device = "cuda", tree_method = "hist", objective = "binary:logistic",
               max_depth = 1, eta = 1), dtrain, nrounds = 1, verbose = 0)
        "cuda"
      }, error = function(e) "cpu")
      assign(".xgb_device", device, envir = .GlobalEnv)
    }
    device <- get(".xgb_device", envir = .GlobalEnv)

    fit <- xgboost::xgb.train(
      params = list(
        device = device,
        tree_method = "hist",
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = params$max_depth,
        eta = params$eta
      ),
      data = dtrain,
      nrounds = params$nrounds,
      verbose = 0
    )

    # Native SHAP via predict(predcontrib = TRUE)
    shap_vals <- predict(fit, dtrain, predcontrib = TRUE)
    # Last column is bias term — remove it
    shap_vals <- shap_vals[, -ncol(shap_vals), drop = FALSE]

    # Mean absolute SHAP as importance
    importance <- colMeans(abs(shap_vals))
    selected <- select_top_k(importance, params$top_k_method, params$top_n)

    list(selected = selected, importance = importance,
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = TRUE, message = paste0("device=", device))
  }, error = function(e) {
    list(selected = rep(FALSE, ncol(X)), importance = rep(NA_real_, ncol(X)),
         time = as.numeric(proc.time()["elapsed"] - t0),
         converged = FALSE, message = e$message)
  })
}

# ==============================================================================
# HELPER: Top-k selection by elbow method
# ==============================================================================

#' Select top features by elbow detection or fixed number
#'
#' @param importance Named numeric vector of importance scores
#' @param method "elbow" or "fixed"
#' @param top_n Fallback number of features if elbow fails
#' @return Logical vector (length p)
select_top_k <- function(importance, method = "elbow", top_n = 20) {
  p <- length(importance)
  if (all(is.na(importance))) return(rep(FALSE, p))

  sorted_idx <- order(importance, decreasing = TRUE)
  sorted_imp <- importance[sorted_idx]

  k <- top_n  # Default

  if (method == "elbow") {
    # Simple elbow: find maximum curvature point
    # Using second derivative of sorted importance curve
    n_cand <- min(p, max(top_n * 3, 50))
    vals <- sorted_imp[seq_len(n_cand)]
    if (length(vals) > 3) {
      # Normalize
      vals_norm <- (vals - min(vals, na.rm = TRUE)) /
        max(vals - min(vals, na.rm = TRUE), 1e-10)
      # Second derivative
      d2 <- diff(diff(vals_norm))
      # Elbow = point of maximum curvature (most negative second derivative)
      elbow_idx <- which.min(d2) + 1
      if (length(elbow_idx) > 0 && elbow_idx > 1) {
        k <- elbow_idx
      }
    }
  }

  k <- min(k, p)
  selected <- rep(FALSE, p)
  selected[sorted_idx[seq_len(k)]] <- TRUE
  selected
}

# ==============================================================================
# DISPATCHER
# ==============================================================================

#' Run a feature selection method by name
#'
#' @param method_name Character, method name from config
#' @param X Numeric matrix (n × p)
#' @param y Binary response (0/1)
#' @param params Named list of method-specific parameters
#' @return Standard FS result list
run_fs_method <- function(method_name, X, y, params) {
  fn <- switch(method_name,
    "wilcoxon_fdr"        = fs_wilcoxon_fdr,
    "fold_change"         = fs_fold_change,
    "volcano"             = fs_volcano,
    "lasso"               = fs_lasso,
    "elastic_net"         = fs_elastic_net,
    "boruta"              = fs_boruta,
    "rf_importance"       = fs_rf_importance,
    "stability_selection" = fs_stability_selection,
    "knockoff"            = fs_knockoff,
    "horseshoe"           = fs_horseshoe,
    "spike_slab"          = fs_spike_slab,
    "shap_xgboost"        = fs_shap_xgboost,
    NULL
  )

  if (is.null(fn)) {
    cli::cli_alert_danger("Unknown method: {method_name}")
    return(list(
      selected = rep(FALSE, ncol(X)),
      importance = rep(NA_real_, ncol(X)),
      time = 0, converged = FALSE,
      message = paste("Unknown method:", method_name)
    ))
  }

  # Ensure X is a matrix with column names
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))

  fn(X, y, params)
}
