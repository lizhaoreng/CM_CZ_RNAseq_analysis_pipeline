#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run orthogroup-level integration of all functional annotation results
#
# Usage:
#   bash 07_run_integration.sh
#
# Requirements:
#   - config.sh in the same directory
#   - 06_integrate_all_annotations.py in the same directory
# ==============================================================================

CONFIG_FILE="./config.sh"
SCRIPT_FILE="./06_integrate_all_annotations.py"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "[ERROR] config.sh not found: ${CONFIG_FILE}" >&2
    exit 1
fi

if [[ ! -f "${SCRIPT_FILE}" ]]; then
    echo "[ERROR] Python integration script not found: ${SCRIPT_FILE}" >&2
    exit 1
fi

source "${CONFIG_FILE}"

mkdir -p "${INTEGRATED_DIR}"
mkdir -p "${LOG_DIR}"

# Export paths for Python script
export WORK_DIR
export OG_CSV
export OG_XLSX
export INTERPRO_DIR
export KOFAM_DIR
export EGGNOG_DIR
export BLAST_DIR
export PHI_DIR
export MEROPS_DIR
export DBCAN_DIR
export INTEGRATED_DIR

echo "============================================================"
echo "Functional annotation integration"
echo "============================================================"
echo "Start time: $(date)"
echo "WORK_DIR: ${WORK_DIR}"
echo "Output directory: ${INTEGRATED_DIR}"
echo "============================================================"

python "${SCRIPT_FILE}" \
    > "${INTEGRATED_DIR}/integration.run.log" \
    2> "${INTEGRATED_DIR}/integration.err.log"

echo "============================================================"
echo "Integration completed"
echo "End time: $(date)"
echo "============================================================"

echo ""
echo "Main output files:"
ls -lh "${INTEGRATED_DIR}" || true

echo ""
echo "Summary statistics:"
if [[ -f "${INTEGRATED_DIR}/integration_summary_statistics.csv" ]]; then
    cat "${INTEGRATED_DIR}/integration_summary_statistics.csv"
else
    echo "[WARNING] integration_summary_statistics.csv was not found."
fi

echo ""
echo "Log files:"
echo " - ${INTEGRATED_DIR}/integration.run.log"
echo " - ${INTEGRATED_DIR}/integration.err.log"
