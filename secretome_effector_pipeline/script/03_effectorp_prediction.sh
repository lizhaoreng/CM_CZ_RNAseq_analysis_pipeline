#!/usr/bin/env bash

# ==============================================================================
# 03_effectorp_prediction.sh
#
# Predict candidate effectors from classical secreted proteins using EffectorP.
# ==============================================================================

set -euo pipefail

CONFIG="${1:-./config.sh}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.sh not found."
    exit 1
fi

source "${CONFIG}"

LOG_FILE="${LOG_DIR}/03_effectorp_prediction.log"

mkdir -p \
    "${SECRETOME_DIR}/${SAMPLE1}/effectorp" \
    "${SECRETOME_DIR}/${SAMPLE2}/effectorp" \
    "${LOG_DIR}"

{
    echo "===== Step 3: EffectorP prediction ====="
    echo "Start time: $(date)"
    echo ""
} > "${LOG_FILE}"

if ! command -v EffectorP >/dev/null 2>&1; then
    echo "ERROR: EffectorP not found in PATH." | tee -a "${LOG_FILE}"
    exit 1
fi

echo "EffectorP path: $(command -v EffectorP)" | tee -a "${LOG_FILE}"

run_effectorp_for_sample() {
    local sample="$1"

    local input_fasta="${SECRETOME_DIR}/${sample}/${sample}_secreted_proteins.faa"
    local out_dir="${SECRETOME_DIR}/${sample}/effectorp"
    local out_tsv="${out_dir}/${sample}_effectorp.tsv"
    local effector_fasta="${out_dir}/${sample}_effectors.faa"
    local noneffector_fasta="${out_dir}/${sample}_noneffectors.faa"

    if [[ ! -s "${input_fasta}" ]]; then
        echo "ERROR: Secreted protein FASTA not found or empty: ${input_fasta}" | tee -a "${LOG_FILE}"
        exit 1
    fi

    echo "" | tee -a "${LOG_FILE}"
    echo "Running EffectorP for sample: ${sample}" | tee -a "${LOG_FILE}"

    EffectorP -f \
        -i "${input_fasta}" \
        -o "${out_tsv}" \
        -E "${effector_fasta}" \
        -N "${noneffector_fasta}"

    local total
    local effectors
    local noneffectors

    total=$(grep -c "^>" "${input_fasta}" || true)
    effectors=$(grep -c "^>" "${effector_fasta}" || true)
    noneffectors=$(grep -c "^>" "${noneffector_fasta}" || true)

    echo "${sample} secreted candidates: ${total}" | tee -a "${LOG_FILE}"
    echo "${sample} predicted effectors: ${effectors}" | tee -a "${LOG_FILE}"
    echo "${sample} predicted non-effectors: ${noneffectors}" | tee -a "${LOG_FILE}"

    echo -e "${sample}\t${total}\t${effectors}\t${noneffectors}"
}

SUMMARY_FILE="${SECRETOME_DIR}/03_effectorp_summary.tsv"

{
    echo -e "Sample\tSecreted_candidates\tEffectors\tNon_effectors"
    run_effectorp_for_sample "${SAMPLE1}"
    run_effectorp_for_sample "${SAMPLE2}"
} > "${SUMMARY_FILE}"

echo "" | tee -a "${LOG_FILE}"
echo "Summary table saved: ${SUMMARY_FILE}" | tee -a "${LOG_FILE}"
cat "${SUMMARY_FILE}" | tee -a "${LOG_FILE}"

echo "Step 3 completed." | tee -a "${LOG_FILE}"
echo "End time: $(date)" | tee -a "${LOG_FILE}"
