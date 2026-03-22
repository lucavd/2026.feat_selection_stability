#!/usr/bin/env bash
# Build manuscript — Pandoc Markdown → DOCX / PDF
# Usage: ./build.sh [docx|pdf|both]  (default: docx)
set -euo pipefail
cd "$(dirname "$0")"

FORMAT="${1:-docx}"
STEM="manuscript"
BIB="references.bib"
CSL="csl/numeric.csl"

COMMON_OPTS=(
  --citeproc
  --bibliography="$BIB"
  --csl="$CSL"
  --number-sections
  --standalone
)

build_docx() {
  echo "Building DOCX..."
  pandoc "$STEM.md" "${COMMON_OPTS[@]}" \
    --reference-doc=reference.docx \
    -o "$STEM.docx" 2>/dev/null || \
  pandoc "$STEM.md" "${COMMON_OPTS[@]}" \
    -o "$STEM.docx"
  echo "  → $STEM.docx"
}

build_pdf() {
  echo "Building PDF..."
  pandoc "$STEM.md" "${COMMON_OPTS[@]}" \
    --pdf-engine=xelatex \
    -V geometry:margin=2.5cm \
    -V fontsize=11pt \
    -V linestretch=1.5 \
    -o "$STEM.pdf"
  echo "  → $STEM.pdf"
}

case "$FORMAT" in
  docx) build_docx ;;
  pdf)  build_pdf ;;
  both) build_docx; build_pdf ;;
  *)    echo "Usage: $0 [docx|pdf|both]"; exit 1 ;;
esac

echo "Done."
