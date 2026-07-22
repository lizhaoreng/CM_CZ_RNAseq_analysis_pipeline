#!/usr/bin/env bash

# ==============================================================================
# Secretome and effector prediction pipeline
#
# Steps:
#   1. Prepare input protein FASTA files and check FASTA format
#   2. Predict signal peptides using SignalP6 and filter proteins without TMHs
#   3. Predict candidate effectors using EffectorP
#   4. Predict subcellular localization of candidate effectors using DeepLoc2
#
# Usage:
#   bash run_secretome_effector_pipeline.sh
# ==============================================================================

set -euo pipefail

CONFIG="./config.sh"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.sh not found."
    exit 1
fi

source "${CONFIG}"

echo "============================================================"
echo "Secretome and effector prediction pipeline"
echo "Start time: $(date)"
echo "============================================================"

echo ""
echo "[Step 1] Preparing protein FASTA input files..."
bash scripts/01_prepare_secretome_input.sh "${CONFIG}"

echo ""
echo "[Step 2] Running SignalP6 and TMHMM filtering..."
bash scripts/02_signalp_tmhmm_filter.sh "${CONFIG}"

echo ""
echo "[Step 3] Running EffectorP prediction..."
bash scripts/03_effectorp_prediction.sh "${CONFIG}"

echo ""
echo "[Step 4] Running DeepLoc2 localization prediction..."
bash scripts/04_deeploc2_effector_localization.sh "${CONFIG}"

echo ""
echo "============================================================"
echo "Pipeline completed successfully"
echo "End time: $(date)"
echo "============================================================"

echo ""
echo "Main output directory:"
echo "  ${SECRETOME_DIR}"
