#!/usr/bin/env bash

# ==============================================================================
# 02_signalp_tmhmm_filter.sh
#
# Predict signal peptides using SignalP6 and remove proteins with predicted
# transmembrane helices using TMHMM.
# ==============================================================================

set -euo pipefail

CONFIG="${1:-./config.sh}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.sh not found."
    exit 1
fi

source "${CONFIG}"

LOG_FILE="${LOG_DIR}/02_signalp_tmhmm_filter.log"

mkdir -p \
    "${SECRETOME_DIR}/${SAMPLE1}/signalp" \
    "${SECRETOME_DIR}/${SAMPLE2}/signalp" \
    "${SECRETOME_DIR}/${SAMPLE1}/tmhmm" \
    "${SECRETOME_DIR}/${SAMPLE2}/tmhmm" \
    "${LOG_DIR}"

{
    echo "===== Step 2: SignalP6 prediction and TMHMM filtering ====="
    echo "Start time: $(date)"
    echo ""
} > "${LOG_FILE}"

SAMPLE1_INPUT="${SECRETOME_DIR}/${SAMPLE1}/${SAMPLE1}_proteins.faa"
SAMPLE2_INPUT="${SECRETOME_DIR}/${SAMPLE2}/${SAMPLE2}_proteins.faa"

for file in "${SAMPLE1_INPUT}" "${SAMPLE2_INPUT}"; do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Input file not found or empty: ${file}" | tee -a "${LOG_FILE}"
        exit 1
    fi
done

if ! command -v signalp6 >/dev/null 2>&1; then
    echo "ERROR: signalp6 not found in PATH." | tee -a "${LOG_FILE}"
    exit 1
fi

if ! command -v tmhmm >/dev/null 2>&1; then
    echo "ERROR: tmhmm not found in PATH." | tee -a "${LOG_FILE}"
    exit 1
fi

echo "signalp6 path: $(command -v signalp6)" | tee -a "${LOG_FILE}"
echo "tmhmm path: $(command -v tmhmm)" | tee -a "${LOG_FILE}"

run_for_sample() {
    local sample="$1"

    local input_fasta="${SECRETOME_DIR}/${sample}/${sample}_proteins.faa"
    local cleaned_fasta="${SECRETOME_DIR}/${sample}/${sample}_proteins_cleaned.faa"
    local signalp_dir="${SECRETOME_DIR}/${sample}/signalp"
    local tmhmm_dir="${SECRETOME_DIR}/${sample}/tmhmm"

    local signalp_result="${signalp_dir}/prediction_results.txt"
    local signalp_fasta="${SECRETOME_DIR}/${sample}/${sample}_signalp_proteins.faa"
    local tmhmm_result="${tmhmm_dir}/tmhmm_results.txt"
    local secreted_fasta="${SECRETOME_DIR}/${sample}/${sample}_secreted_proteins.faa"

    echo "" | tee -a "${LOG_FILE}"
    echo "Processing sample: ${sample}" | tee -a "${LOG_FILE}"

    echo "Cleaning protein FASTA..." | tee -a "${LOG_FILE}"
    python3 scripts/clean_protein_fasta.py \
        --input "${input_fasta}" \
        --output "${cleaned_fasta}" | tee -a "${LOG_FILE}"

    local total_count
    local cleaned_count

    total_count=$(grep -c "^>" "${input_fasta}" || true)
    cleaned_count=$(grep -c "^>" "${cleaned_fasta}" || true)

    echo "Total proteins: ${total_count}" | tee -a "${LOG_FILE}"
    echo "Cleaned proteins: ${cleaned_count}" | tee -a "${LOG_FILE}"

    echo "Running SignalP6..." | tee -a "${LOG_FILE}"
    rm -rf "${signalp_dir:?}/"*
    signalp6 \
        --fastafile "${cleaned_fasta}" \
        --organism "${SIGNALP_ORGANISM}" \
        --format "${SIGNALP_FORMAT}" \
        --output_dir "${signalp_dir}"

    if [[ ! -s "${signalp_result}" ]]; then
        echo "ERROR: SignalP6 result not found or empty: ${signalp_result}" | tee -a "${LOG_FILE}"
        exit 1
    fi

    echo "Extracting SignalP-positive proteins..." | tee -a "${LOG_FILE}"
    python3 scripts/extract_secretome_candidates.py \
        --result "${signalp_result}" \
        --fasta "${cleaned_fasta}" \
        --output "${signalp_fasta}" \
        --mode signalp | tee -a "${LOG_FILE}"

    local signalp_count
    signalp_count=$(grep -c "^>" "${signalp_fasta}" || true)

    echo "SignalP-positive proteins: ${signalp_count}" | tee -a "${LOG_FILE}"

    if [[ "${signalp_count}" -eq 0 ]]; then
        echo "ERROR: No SignalP-positive proteins detected for ${sample}." | tee -a "${LOG_FILE}"
        exit 1
    fi

    echo "Running TMHMM..." | tee -a "${LOG_FILE}"
    tmhmm "${signalp_fasta}" > "${tmhmm_result}"

    if [[ ! -s "${tmhmm_result}" ]]; then
        echo "ERROR: TMHMM result not found or empty: ${tmhmm_result}" | tee -a "${LOG_FILE}"
        exit 1
    fi

    echo "Extracting proteins without predicted transmembrane helices..." | tee -a "${LOG_FILE}"
    python3 scripts/extract_secretome_candidates.py \
        --result "${tmhmm_result}" \
        --fasta "${signalp_fasta}" \
        --output "${secreted_fasta}" \
        --mode tmhmm | tee -a "${LOG_FILE}"

    local secreted_count
    secreted_count=$(grep -c "^>" "${secreted_fasta}" || true)

    echo "Secreted candidates: ${secreted_count}" | tee -a "${LOG_FILE}"

    echo -e "${sample}\t${total_count}\t${cleaned_count}\t${signalp_count}\t${secreted_count}"
}

SUMMARY_FILE="${SECRETOME_DIR}/02_signalp_tmhmm_summary.tsv"

{
    echo -e "Sample\tTotal_proteins\tCleaned_proteins\tSignalP_positive\tSecreted_candidates"
    run_for_sample "${SAMPLE1}"
    run_for_sample "${SAMPLE2}"
} > "${SUMMARY_FILE}"

echo "" | tee -a "${LOG_FILE}"
echo "Summary table saved: ${SUMMARY_FILE}" | tee -a "${LOG_FILE}"
cat "${SUMMARY_FILE}" | tee -a "${LOG_FILE}"

echo "Step 2 completed." | tee -a "${LOG_FILE}"
echo "End time: $(date)" | tee -a "${LOG_FILE}"
