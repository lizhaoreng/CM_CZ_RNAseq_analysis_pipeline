#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run all fungal functional annotation steps
# ==============================================================================

source ./config.sh

echo "[$(date '+%F %T')] Functional annotation pipeline started"

bash 01_run_interpro_kofam_eggnog.sh -g both -s all -t "${THREADS}"
bash 02_run_ncbi_fungi_blast.sh
bash 03_run_dbcan_cwde.sh
bash 04_run_merops.sh
bash 05_run_phi_uniprot.sh
bash 07_run_integration.sh

export WORK_DIR
python 08_build_cwde_4class_table.py \
  > "${DBCAN_DIR}/cwde_4class_table.run.log" \
  2> "${DBCAN_DIR}/cwde_4class_table.err.log"

echo "[$(date '+%F %T')] Functional annotation pipeline completed"
