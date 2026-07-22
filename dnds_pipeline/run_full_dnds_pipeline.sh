#!/usr/bin/env bash

# ==============================================================================
# Full dN/dS analysis pipeline
#
# Steps:
#   1. Prepare and rename orthologous CDS/protein sequences
#   2. Extract one-to-one orthologous sequence pairs
#   3. Run protein alignment, codon back-translation, and codeml
#   4. Calculate alignment quality statistics
#   5. Parse codeml results and apply quality control
#   6. Analyze evolutionary rates of effectors, CWDEs, and background genes
#
# Usage:
#   bash run_full_dnds_pipeline.sh
# ==============================================================================

set -euo pipefail

CONFIG="config.yaml"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.yaml not found."
    exit 1
fi

echo "============================================================"
echo "Full dN/dS evolutionary analysis pipeline"
echo "Start time: $(date)"
echo "============================================================"

echo ""
echo "[Step 1] Preparing ortholog data and renaming FASTA sequences..."
python3 scripts/01_prepare_ortholog_data.py --config "${CONFIG}"

echo ""
echo "[Step 2] Extracting one-to-one orthologous sequence pairs..."
python3 scripts/02_extract_ortholog_pairs.py --config "${CONFIG}"

echo ""
echo "[Step 3] Running codon alignment and codeml pairwise dN/dS analysis..."
bash scripts/03_run_codon_alignment_codeml.sh "${CONFIG}"

echo ""
echo "[Step 3b] Calculating codon alignment quality statistics..."
python3 scripts/03b_alignment_quality_stats.py --config "${CONFIG}"

echo ""
echo "[Step 4] Parsing codeml pairwise dN/dS results..."
python3 scripts/04_parse_codeml_results.py --config "${CONFIG}"

echo ""
echo "[Step 5] Analyzing evolutionary rates by functional gene category..."
Rscript scripts/05_analyze_dnds_functional_categories.R "${CONFIG}"

echo ""
echo "============================================================"
echo "Pipeline completed successfully"
echo "End time: $(date)"
echo "============================================================"

echo ""
echo "Main outputs:"
echo "  - results/dnds_enhanced_complete.csv"
echo "  - results/dnds_enhanced_clean.csv"
echo "  - results/alignment_length_stats.csv"
echo "  - results/evolutionary_analysis/"
