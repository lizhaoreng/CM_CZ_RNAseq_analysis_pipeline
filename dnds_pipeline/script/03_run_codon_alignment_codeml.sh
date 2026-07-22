#!/usr/bin/env bash

# ==============================================================================
# 03_run_codon_alignment_codeml.sh
#
# Run protein alignment, codon back-translation, and pairwise codeml analysis.
#
# Requirements:
#   mafft
#   pal2nal.pl
#   codeml
# ==============================================================================

set -euo pipefail

CONFIG="${1:-config.yaml}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.yaml not found."
    exit 1
fi

get_yaml_value() {
    python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c$1)"
}

SPECIES1_CODE=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['species']['species1_code'])")
SPECIES2_CODE=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['species']['species2_code'])")

PAIR_PROT_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['pair_protein_dir'])")
PAIR_CDS_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['pair_cds_dir'])")
CODON_ALN_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['codon_alignment_dir'])")
CODEML_OUT_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['codeml_out_dir'])")
LOG_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['logs_dir'])")

CODON_FREQ=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['codeml']['codon_freq'])")
OMEGA_INITIAL=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['codeml']['omega_initial'])")
KAPPA_INITIAL=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['codeml']['kappa_initial'])")
CLEANDATA=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['codeml']['cleandata'])")

TEMP_DIR="results/temp_codeml"

mkdir -p "${CODON_ALN_DIR}" "${CODEML_OUT_DIR}" "${LOG_DIR}" "${TEMP_DIR}"

if [[ ! -d "${PAIR_PROT_DIR}" || ! -d "${PAIR_CDS_DIR}" ]]; then
    echo "ERROR: Missing input directories: ${PAIR_PROT_DIR} or ${PAIR_CDS_DIR}"
    exit 1
fi

PROT_COUNT=$(find "${PAIR_PROT_DIR}" -name "*.faa" | wc -l)

echo "============================================================"
echo "Step 3: Codon alignment and pairwise codeml analysis"
echo "============================================================"
echo "Species: ${SPECIES1_CODE} vs ${SPECIES2_CODE}"
echo "Protein pair files: ${PROT_COUNT}"

if [[ "${PROT_COUNT}" -eq 0 ]]; then
    echo "ERROR: No protein pair files found."
    exit 1
fi

SUCCESS_COUNT=0
FAILED_ALIGNMENT=0
FAILED_BACKTRANSLATION=0
FAILED_CODEML=0
TOTAL_COUNT=0

: > "${LOG_DIR}/codeml_success.log"
: > "${LOG_DIR}/codeml_failed.log"

for PROT_FILE in "${PAIR_PROT_DIR}"/*.faa; do
    OG=$(basename "${PROT_FILE}" .faa)
    CDS_FILE="${PAIR_CDS_DIR}/${OG}.cds"

    TOTAL_COUNT=$((TOTAL_COUNT + 1))

    if [[ ! -f "${CDS_FILE}" ]]; then
        echo "FAIL:${OG}:NO_CDS" >> "${LOG_DIR}/codeml_failed.log"
        continue
    fi

    PROT_ALN="${TEMP_DIR}/${OG}.protein.aln"
    CODON_ALN="${CODON_ALN_DIR}/${OG}.phy"
    TREE_FILE="${TEMP_DIR}/${OG}.tree"
    CTL_FILE="${TEMP_DIR}/${OG}.ctl"
    OUT_FILE="${CODEML_OUT_DIR}/${OG}.out"

    if ! mafft --auto --quiet "${PROT_FILE}" > "${PROT_ALN}" 2>/dev/null; then
        FAILED_ALIGNMENT=$((FAILED_ALIGNMENT + 1))
        echo "FAIL:${OG}:MAFFT" >> "${LOG_DIR}/codeml_failed.log"
        continue
    fi

    if ! pal2nal.pl "${PROT_ALN}" "${CDS_FILE}" -output paml > "${CODON_ALN}" 2>/dev/null; then
        FAILED_BACKTRANSLATION=$((FAILED_BACKTRANSLATION + 1))
        echo "FAIL:${OG}:PAL2NAL" >> "${LOG_DIR}/codeml_failed.log"
        rm -f "${PROT_ALN}"
        continue
    fi

    echo "(${SPECIES1_CODE}_${OG}, ${SPECIES2_CODE}_${OG});" > "${TREE_FILE}"

    cat > "${CTL_FILE}" << EOF
      seqfile = ${CODON_ALN}
     treefile = ${TREE_FILE}
      outfile = ${OUT_FILE}
        noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = ${CODON_FREQ}
        model = 0
      NSsites = 0
        icode = 0
    fix_kappa = 0
        kappa = ${KAPPA_INITIAL}
    fix_omega = 0
        omega = ${OMEGA_INITIAL}
    cleandata = ${CLEANDATA}
EOF

    codeml "${CTL_FILE}" >/dev/null 2>&1 || true

    if [[ -s "${OUT_FILE}" ]] && grep -qi "omega" "${OUT_FILE}"; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        echo "SUCCESS:${OG}" >> "${LOG_DIR}/codeml_success.log"
    else
        FAILED_CODEML=$((FAILED_CODEML + 1))
        echo "FAIL:${OG}:CODEML" >> "${LOG_DIR}/codeml_failed.log"
    fi

    rm -f "${PROT_ALN}" "${TREE_FILE}" "${CTL_FILE}"

    if [[ $((TOTAL_COUNT % 500)) -eq 0 ]]; then
        echo "Processed ${TOTAL_COUNT}/${PROT_COUNT}; success: ${SUCCESS_COUNT}"
    fi
done

rm -rf "${TEMP_DIR}"

echo ""
echo "Step 3 completed."
echo "Total processed: ${TOTAL_COUNT}"
echo "Successful codeml runs: ${SUCCESS_COUNT}"
echo "MAFFT failures: ${FAILED_ALIGNMENT}"
echo "PAL2NAL failures: ${FAILED_BACKTRANSLATION}"
echo "codeml failures: ${FAILED_CODEML}"

if [[ "${SUCCESS_COUNT}" -eq 0 ]]; then
    echo "ERROR: No successful codeml results."
    exit 1
fi
