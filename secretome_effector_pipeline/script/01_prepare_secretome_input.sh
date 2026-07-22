#!/usr/bin/env bash

# ==============================================================================
# 01_prepare_secretome_input.sh
#
# Prepare protein FASTA files and check FASTA format.
# ==============================================================================

set -euo pipefail

CONFIG="${1:-./config.sh}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.sh not found."
    exit 1
fi

source "${CONFIG}"

mkdir -p \
    "${SECRETOME_DIR}/${SAMPLE1}" \
    "${SECRETOME_DIR}/${SAMPLE2}" \
    "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/01_prepare_secretome_input.log"

{
    echo "===== Step 1: Prepare protein FASTA input files ====="
    echo "Start time: $(date)"
    echo "Project directory: ${PROJECT_DIR}"
    echo "Reference directory: ${REF_DIR}"
    echo ""
} > "${LOG_FILE}"

echo "Checking input protein FASTA files..." | tee -a "${LOG_FILE}"

if [[ ! -f "${SAMPLE1_PROTEINS}" ]]; then
    echo "ERROR: Protein FASTA file not found: ${SAMPLE1_PROTEINS}" | tee -a "${LOG_FILE}"
    exit 1
fi

if [[ ! -f "${SAMPLE2_PROTEINS}" ]]; then
    echo "ERROR: Protein FASTA file not found: ${SAMPLE2_PROTEINS}" | tee -a "${LOG_FILE}"
    exit 1
fi

echo "Found ${SAMPLE1} protein FASTA: ${SAMPLE1_PROTEINS}" | tee -a "${LOG_FILE}"
echo "Found ${SAMPLE2} protein FASTA: ${SAMPLE2_PROTEINS}" | tee -a "${LOG_FILE}"

cp "${SAMPLE1_PROTEINS}" "${SECRETOME_DIR}/${SAMPLE1}/${SAMPLE1}_proteins.faa"
cp "${SAMPLE2_PROTEINS}" "${SECRETOME_DIR}/${SAMPLE2}/${SAMPLE2}_proteins.faa"

SAMPLE1_COUNT=$(grep -c "^>" "${SECRETOME_DIR}/${SAMPLE1}/${SAMPLE1}_proteins.faa" || true)
SAMPLE2_COUNT=$(grep -c "^>" "${SECRETOME_DIR}/${SAMPLE2}/${SAMPLE2}_proteins.faa" || true)

echo "${SAMPLE1} protein count: ${SAMPLE1_COUNT}" | tee -a "${LOG_FILE}"
echo "${SAMPLE2} protein count: ${SAMPLE2_COUNT}" | tee -a "${LOG_FILE}"

echo "Checking FASTA format..." | tee -a "${LOG_FILE}"

python3 scripts/check_protein_format.py \
    --input "${SECRETOME_DIR}/${SAMPLE1}/${SAMPLE1}_proteins.faa" \
    --log "${LOG_DIR}/${SAMPLE1}_format_check.log"

python3 scripts/check_protein_format.py \
    --input "${SECRETOME_DIR}/${SAMPLE2}/${SAMPLE2}_proteins.faa" \
    --log "${LOG_DIR}/${SAMPLE2}_format_check.log"

README_FILE="${SECRETOME_DIR}/README.md"

cat > "${README_FILE}" << EOF
# Secretome and effector prediction results

This directory contains outputs from the secretome and effector prediction pipeline.

## Samples

- ${SAMPLE1}
- ${SAMPLE2}

## Main steps

1. Protein FASTA format checking
2. SignalP6 signal peptide prediction
3. TMHMM transmembrane helix filtering
4. EffectorP candidate effector prediction
5. DeepLoc2 localization prediction
EOF

echo "Step 1 completed." | tee -a "${LOG_FILE}"
echo "Output directory: ${SECRETOME_DIR}" | tee -a "${LOG_FILE}"
echo "End time: $(date)" | tee -a "${LOG_FILE}"
