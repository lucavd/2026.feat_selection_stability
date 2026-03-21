#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

mkdir -p logs

echo "[1/2] Starting Telegram watcher..."
Rscript R/09_telegram_progress.R --watch > logs/09_telegram_watch.log 2>&1 &
WATCH_PID=$!
echo "Watcher PID: $WATCH_PID"

echo "[2/2] Starting feature selection (Script 05)..."
# CUDA_VISIBLE_DEVICES="" MUST be set before R starts, not inside R.
# mclapply forks inherit parent's CUDA file descriptors, causing driver
# deadlocks when 40+ workers compete for /dev/nvidia0. Setting this env
# var before launch prevents xgboost from initializing CUDA at all.
CUDA_VISIBLE_DEVICES="" Rscript R/05_feature_selection.R | tee logs/05_feature_selection.log

kill "$WATCH_PID" >/dev/null 2>&1 || true

echo "Feature selection finished."
