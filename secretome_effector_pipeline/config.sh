#!/usr/bin/env bash

# ==============================================================================
# Configuration file for secretome and effector prediction pipeline
# ==============================================================================

# Sample names
SAMPLE1="CM"
SAMPLE2="CZ"

# Project root directory
PROJECT_DIR="$(pwd)"

# Input protein directory
REF_DIR="${PROJECT_DIR}/ref"

# Input protein FASTA files
SAMPLE1_PROTEINS="${REF_DIR}/CM_proteins_cleaned.faa"
SAMPLE2_PROTEINS="${REF_DIR}/CZ_proteins_cleaned.faa"

# Output directories
RESULT_DIR="${PROJECT_DIR}/results"
SECRETOME_DIR="${RESULT_DIR}/secretome"
LOG_DIR="${SECRETOME_DIR}/logs"

# SignalP parameters
SIGNALP_ORGANISM="euk"
SIGNALP_FORMAT="txt"

# DeepLoc2 parameters
DEEPLOC2_MODE="Fast"
DEEPLOC2_DEVICE="cpu"
