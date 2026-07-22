#!/usr/bin/env bash

# ==============================================================================
# Configuration file for fungal protein functional annotation pipeline
# ==============================================================================

# Project directory
WORK_DIR="/path/to/project"

# Input protein directory
INPUT_DIR="${WORK_DIR}/orthofinder_input"

# Protein FASTA files
CM_PROTEINS="${INPUT_DIR}/CM_proteins_cleaned.faa"
CZ_PROTEINS="${INPUT_DIR}/CZ_proteins_cleaned.faa"

# Orthogroup table
OG_CSV="${WORK_DIR}/CM_CZ_OG.expanded.csv"
OG_XLSX="${WORK_DIR}/Complete_CM_CZ_Ortholog_Analysis.xlsx"

# Output directories
ANNOTATION_DIR="${WORK_DIR}/annotation_results"
INTERPRO_DIR="${ANNOTATION_DIR}/interpro"
KOFAM_DIR="${ANNOTATION_DIR}/kofam"
EGGNOG_DIR="${ANNOTATION_DIR}/eggnog"
BLAST_DIR="${INPUT_DIR}/blast_results"
DBCAN_DIR="${WORK_DIR}/dbcan_annotation"
MEROPS_DIR="${WORK_DIR}/merops_annotation"
PHI_DIR="${WORK_DIR}/phi_annotation"
INTEGRATED_DIR="${WORK_DIR}/integrated_annotation_results"
LOG_DIR="${WORK_DIR}/logs"

# Databases
KOFAM_DB="/path/to/kofam_db"
EGGNOG_DATA_DIR="/path/to/eggnog_db"
INTERPROSCAN="/path/to/interproscan.sh"

NCBI_FUNGI_DMND="/path/to/refseq_fungi_complete.dmnd"

DBCAN_DB_DIR="/path/to/dbcan_db"
DBCAN_HMM_DB="${DBCAN_DB_DIR}/dbCAN-HMMdb-V13.txt"
CAZY_DB_FA="${DBCAN_DB_DIR}/CAZyDB.07142024.fa"
CAZY_DB_DMND="${DBCAN_DB_DIR}/CAZyDB.07142024.fa.dmnd"

MEROPS_DMND="/path/to/MEROPS/pepunit_clean.lib.dmnd"

PHI_BASE_DMND="/path/to/phi-base.dmnd"
UNIPROT_FUNGI_DMND="/path/to/uniprot_fungi.dmnd"

# Temporary directory
TMP_BASE="${WORK_DIR}/tmp"

# General parameters
THREADS=64
THREADS_SMALL=16

# DIAMOND parameters
EVALUE=1e-5
MAX_TARGET_SEQS_NCBI=10
MAX_TARGET_SEQS_MEROPS=5
MAX_TARGET_SEQS_PHI=5
MAX_TARGET_SEQS_UNIPROT=3

# Filtering thresholds
MIN_IDENTITY_DBCAN=30
MIN_IDENTITY_MEROPS=30
MIN_IDENTITY_PHI=40
MIN_IDENTITY_UNIPROT=50

# Create directories
mkdir -p \
  "${INTERPRO_DIR}" \
  "${KOFAM_DIR}" \
  "${EGGNOG_DIR}" \
  "${BLAST_DIR}" \
  "${DBCAN_DIR}" \
  "${MEROPS_DIR}" \
  "${PHI_DIR}" \
  "${INTEGRATED_DIR}" \
  "${LOG_DIR}" \
  "${TMP_BASE}"
