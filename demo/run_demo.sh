#!/usr/bin/env bash
#
# run_demo.sh - end-to-end demo runner for the submission package. 
# Runs every analysis step against the small synthetic dataset shipped under demo/data/.
#
# Usage:
#   bash run_demo.sh                       # full demo including CellChat
#   RUN_CELLCHAT=0 bash run_demo.sh        # skip CellChat step
#   R_BIN=... bash run_demo.sh
#                                            # override R interpreters
#
# Expected total runtime on a "normal" desktop: ~5 minutes
#   - demo.Seurat.R            ~90 s
#   - demo.scRNA_pipeline.R    ~60 s
#   - demo.edgeR.R             <5 s
#   - demo.cellchat.R          ~2-3 min (LR inference + 5 pathway plots)
#
# All outputs land under demo/results/ with one subdirectory per step.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DEMO_DIR="$ROOT_DIR/demo"
RESULTS_DIR="$DEMO_DIR/results"
mkdir -p "$RESULTS_DIR"

R_BIN="${R_BIN:-Rscript}"
if ! command -v "$R_BIN" >/dev/null 2>&1; then
  echo "ERROR: $R_BIN not found. Set R_BIN env var to your Rscript path." >&2
  exit 1
fi

RUN_CELLCHAT="${RUN_CELLCHAT:-auto}"

MERGED_INPUT="$DEMO_DIR/data/merged"
RAW_INPUT="$DEMO_DIR/data/raw"
COUNT_FILE="$DEMO_DIR/data/Myeloid.exp.xls"

if [ ! -d "$MERGED_INPUT/Control" ] || [ ! -d "$MERGED_INPUT/PEDV" ]; then
  echo "ERROR: missing merged 10X data under $MERGED_INPUT" >&2
  echo "       re-run gen_demo_data.py + merge_10x.py to recreate it." >&2
  exit 1
fi
if [ ! -f "$COUNT_FILE" ]; then
  echo "ERROR: missing $COUNT_FILE" >&2
  exit 1
fi

echo "=== demo.Seurat.R (CCA integration, scripts 1.Seurat.R equivalent) ==="
SEURAT_OUT="$RESULTS_DIR/seurat"
mkdir -p "$SEURAT_OUT"
"$R_BIN" "$DEMO_DIR/code/demo.Seurat.R" "$MERGED_INPUT" "$SEURAT_OUT"

echo
echo "=== demo.scRNA_pipeline.R (harmony integration, scripts 4.scRNA_pipeline.R equivalent) ==="
PIPE_OUT="$RESULTS_DIR/scRNA_pipeline"
mkdir -p "$PIPE_OUT"
"$R_BIN" "$DEMO_DIR/code/demo.scRNA_pipeline.R" "$DEMO_DIR/config/demo_pipeline_config.R"

echo
echo "=== demo.edgeR.R (scripts 3.edgeR.R equivalent) ==="
EDGER_OUT="$RESULTS_DIR/edger"
mkdir -p "$EDGER_OUT"
cd "$EDGER_OUT"
"$R_BIN" "$DEMO_DIR/code/demo.edgeR.R" --count-file "$COUNT_FILE" --out-prefix "./Myeloid.demo"
cd "$SCRIPT_DIR"

# ---- CellChat (optional, separate R env) ----
cellchat_can_run="yes"
if [ "$RUN_CELLCHAT" = "0" ]; then
  cellchat_can_run="no"
  echo
  echo "=== CellChat step skipped (RUN_CELLCHAT=0) ==="
elif [ "$RUN_CELLCHAT" = "auto" ] && ! command -v "$R_BIN" >/dev/null 2>&1; then
  cellchat_can_run="no"
  echo
  echo "=== CellChat step skipped ($R_BIN not found; set R_BIN to enable) ==="
elif ! command -v "$R_BIN" >/dev/null 2>&1; then  
  cellchat_can_run="no"
  echo
  echo "=== CellChat step skipped ($R_BIN not found) ==="
fi

if [ "$cellchat_can_run" = "yes" ]; then
  echo
  echo "=== demo.cellchat_init.R (CellChat init from Seurat, R=$R_BIN) ==="
  CC_OUT="$RESULTS_DIR/cellchat"
  mkdir -p "$CC_OUT"
  if [ ! -f "$SEURAT_OUT/obj.Rda" ]; then
    echo "WARNING: $SEURAT_OUT/obj.Rda missing; CellChat step aborted." >&2
  else
    "$R_BIN" "$DEMO_DIR/code/demo.cellchat_init.R" "$SEURAT_OUT/obj.Rda" "$CC_OUT"
    echo
    echo "=== demo.cellchat_sample.R (CellChat inference, R=$R_BIN) ==="
    "$R_BIN" "$DEMO_DIR/code/demo.cellchat_sample.R" "$CC_OUT"
  fi
fi

echo
echo "=========================================="
echo "Demo complete."
echo "Inspect:"
echo "  - $SEURAT_OUT/UMAP.pdf, TSNE.pdf, run_summary.txt"
echo "  - $PIPE_OUT/demo_result.seurat.rds, demo_result.UMAP.pdf"
echo "  - $EDGER_OUT/Myeloid.demo.summary.txt"
if [ "$cellchat_can_run" = "yes" ]; then
  echo "  - $CC_OUT/cellchat_summary.txt, CellChat.LR.*.xls, Pathway.*.pdf"
fi
echo "=========================================="
