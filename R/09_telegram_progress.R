# ==============================================================================
# Script 09 — Telegram progress notifier for pipeline
# ==============================================================================
# Purpose: Send progress updates to Telegram, mainly for Script 05 outputs.
# Usage:
#   Rscript R/09_telegram_progress.R             # one snapshot message
#   Rscript R/09_telegram_progress.R --watch     # periodic updates until stopped
#
# Required .env variables at project root:
#   TELEGRAM_BOT_TOKEN=...
#   TELEGRAM_CHAT_ID=...
# Optional:
#   TELEGRAM_INTERVAL_SEC=300
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(cli)
  library(telegram.bot)
})

load_env <- function(env_path = here::here(".env")) {
  if (!file.exists(env_path)) {
    cli::cli_abort("Missing .env file at: {env_path}")
  }
  readRenviron(env_path)
}

build_progress <- function(config, start_time, start_done) {
  sim_dir <- here::here(config$paths$simulated_data)
  fs_dir <- here::here(config$paths$results_fs)

  n_sim <- length(list.files(
    sim_dir, pattern = "^sim_.*\\.rds$", full.names = FALSE
  ))
  n_methods <- length(config$methods)
  n_total <- n_sim * n_methods
  n_done <- length(list.files(
    fs_dir, pattern = "^fs_.*\\.rds$", full.names = FALSE
  ))

  pct <- if (n_total > 0) 100 * n_done / n_total else NA_real_

  # ETA based on rate since watcher started
  elapsed_min <- as.numeric(
    difftime(Sys.time(), start_time, units = "mins")
  )
  new_done <- n_done - start_done
  rate <- if (elapsed_min > 1 && new_done > 0) {
    new_done / elapsed_min
  } else NA_real_
  eta_min <- if (!is.na(rate) && rate > 0) {
    round((n_total - n_done) / rate)
  } else NA_integer_

  list(
    n_sim = n_sim,
    n_methods = n_methods,
    n_total = n_total,
    n_done = n_done,
    n_left = max(n_total - n_done, 0L),
    pct = pct,
    rate = rate,
    eta_min = eta_min,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

format_message <- function(prog) {
  pct_txt <- if (is.na(prog$pct)) {
    "NA"
  } else {
    sprintf("%.1f", prog$pct)
  }
  rate_txt <- if (is.na(prog$rate)) {
    "calculating..."
  } else {
    sprintf("%.1f jobs/min", prog$rate)
  }
  eta_txt <- if (is.na(prog$eta_min)) {
    "calculating..."
  } else if (prog$eta_min > 60) {
    sprintf("%dh %dm", prog$eta_min %/% 60, prog$eta_min %% 60)
  } else {
    sprintf("%d min", prog$eta_min)
  }

  paste0(
    "*Feature Selection Progress*\n",
    "Time: `", prog$timestamp, "`\n",
    "Completed: `", prog$n_done, " / ", prog$n_total,
    "` (", pct_txt, "%)\n",
    "Rate: `", rate_txt, "`\n",
    "ETA: `", eta_txt, "`\n",
    "Remaining: `", prog$n_left, "`"
  )
}

send_message <- function(bot, chat_id, text) {
  tryCatch(
    {
      bot$sendMessage(chat_id = chat_id, text = text, parse_mode = "Markdown")
      TRUE
    },
    error = function(e) {
      cli::cli_alert_warning("Telegram send failed: {e$message}")
      FALSE
    }
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  watch_mode <- "--watch" %in% args

  load_env()

  token <- Sys.getenv("TELEGRAM_BOT_TOKEN", unset = "")
  chat_id <- Sys.getenv("TELEGRAM_CHAT_ID", unset = "")
  interval_sec <- suppressWarnings(as.integer(Sys.getenv("TELEGRAM_INTERVAL_SEC", unset = "300")))
  if (is.na(interval_sec) || interval_sec < 30L) interval_sec <- 300L

  if (token == "" || chat_id == "") {
    cli::cli_abort("Set TELEGRAM_BOT_TOKEN and TELEGRAM_CHAT_ID in .env")
  }

  cfg_path <- here::here("config", "config.yaml")
  if (!file.exists(cfg_path)) {
    cli::cli_abort("Missing config file: {cfg_path}")
  }
  config <- yaml::read_yaml(cfg_path)

  bot <- telegram.bot::Bot(token = token)

  start_time <- Sys.time()
  fs_dir <- here::here(config$paths$results_fs)
  start_done <- length(list.files(
    fs_dir, pattern = "^fs_.*\\.rds$"
  ))

  if (!watch_mode) {
    prog <- build_progress(config, start_time, start_done)
    ok <- send_message(bot, chat_id, format_message(prog))
    if (ok) cli::cli_alert_success("Snapshot sent")
    return(invisible(NULL))
  }

  cli::cli_alert_info(
    "Watch mode ON (every {interval_sec}s). Ctrl+C to stop"
  )
  send_message(
    bot, chat_id,
    "*Telegram watcher started*\nPolling every {interval_sec}s"
  )

  last_done <- -1L
  repeat {
    prog <- build_progress(config, start_time, start_done)

    if (prog$n_done != last_done) {
      send_message(bot, chat_id, format_message(prog))
      last_done <- prog$n_done
    }

    Sys.sleep(interval_sec)
  }
}

main()
