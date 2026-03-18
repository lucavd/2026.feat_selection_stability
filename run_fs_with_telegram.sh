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
Rscript R/05_feature_selection.R | tee logs/05_feature_selection.log

kill "$WATCH_PID" >/dev/null 2>&1 || true

echo "Feature selection finished."
