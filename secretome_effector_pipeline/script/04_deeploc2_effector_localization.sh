#!/usr/bin/env bash

# ==============================================================================
# 04_deeploc2_effector_localization.sh
#
# Predict subcellular localization of EffectorP-predicted candidate effectors
# using DeepLoc2.
# ==============================================================================

set -euo pipefail

CONFIG="${1:-./config.sh}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.sh not found."
    exit 1
fi

source "${CONFIG}"

LOG_FILE="${LOG_DIR}/04_deeploc2_effector_localization.log"

mkdir -p \
    "${SECRETOME_DIR}/${SAMPLE1}/deeploc2_effector" \
    "${SECRETOME_DIR}/${SAMPLE2}/deeploc2_effector" \
    "${LOG_DIR}"

{
    echo "===== Step 4: DeepLoc2 localization prediction for candidate effectors ====="
    echo "Start time: $(date)"
    echo ""
} > "${LOG_FILE}"

if ! command -v deeploc2 >/dev/null 2>&1; then
    echo "ERROR: deeploc2 not found in PATH." | tee -a "${LOG_FILE}"
    exit 1
fi

echo "deeploc2 path: $(command -v deeploc2)" | tee -a "${LOG_FILE}"

run_deeploc2_for_sample() {
    local sample="$1"

    local input_fasta="${SECRETOME_DIR}/${sample}/effectorp/${sample}_effectors.faa"
    local output_dir="${SECRETOME_DIR}/${sample}/deeploc2_effector"

    if [[ ! -s "${input_fasta}" ]]; then
        echo "ERROR: Effector FASTA not found or empty: ${input_fasta}" | tee -a "${LOG_FILE}"
        exit 1
    fi

    echo "" | tee -a "${LOG_FILE}"
    echo "Running DeepLoc2 for sample: ${sample}" | tee -a "${LOG_FILE}"

    deeploc2 \
        -f "${input_fasta}" \
        -o "${output_dir}" \
        -m "${DEEPLOC2_MODE}" \
        -d "${DEEPLOC2_DEVICE}"

    echo "${sample} DeepLoc2 output: ${output_dir}" | tee -a "${LOG_FILE}"
}

run_deeploc2_for_sample "${SAMPLE1}"
run_deeploc2_for_sample "${SAMPLE2}"

echo "" | tee -a "${LOG_FILE}"
echo "Step 4 completed." | tee -a "${LOG_FILE}"
echo "End time: $(date)" | tee -a "${LOG_FILE}"
