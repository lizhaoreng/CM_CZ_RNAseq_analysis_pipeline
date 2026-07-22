#!/usr/bin/env bash

# ==============================================================================
# Configuration file for RNA-seq preprocessing, mapping, splitting, and counting
# ==============================================================================

# Number of threads
THREADS_FASTQC=4
THREADS_TRIMMOMATIC=8
THREADS_HISAT2_BUILD=20
THREADS_HISAT2=20
THREADS_SAMTOOLS=10
THREADS_FEATURECOUNTS=20

# Maximum number of parallel FastQC jobs
MAX_FASTQC_JOBS=24

# Project root directory
PROJECT_DIR="/path/to/project"

# Raw data directory
RAW_DATA_DIR="${PROJECT_DIR}/Data"

# Output directories
FASTQC_PRETRIM_DIR="${PROJECT_DIR}/fastqc_pretrim"
TRIMMED_DIR="${PROJECT_DIR}/trimmed"
FASTQC_POSTTRIM_DIR="${PROJECT_DIR}/fastqc_posttrim"
GENOME_DIR="${PROJECT_DIR}/genome"
INDEX_DIR="${GENOME_DIR}/index"
MAPPING_DIR="${PROJECT_DIR}/mapping"
BAM_DIR="${MAPPING_DIR}/bam"
MAPPING_STAT_DIR="${MAPPING_DIR}/stat"
SPLIT_BAM_DIR="${MAPPING_DIR}/split_bam"
QUANT_DIR="${PROJECT_DIR}/featureCounts"
LOG_DIR="${PROJECT_DIR}/logs"

# Reference genome FASTA files
MAIZE_FA="${GENOME_DIR}/Zea_mays_B73v5.fa"
CM_FA="${GENOME_DIR}/CM_raw.fa"
CZ_FA="${GENOME_DIR}/CZ_raw.fa"

# Reference annotation GFF3 files
MAIZE_GFF3="${GENOME_DIR}/ZM_annotation.gff3"
CM_GFF3="${GENOME_DIR}/CM_annotation.gff3"
CZ_GFF3="${GENOME_DIR}/CZ_annotation.gff3"

# Converted GTF files
MAIZE_GTF="${GENOME_DIR}/ZM_annotation.gtf"
CM_GTF="${GENOME_DIR}/CM_annotation.gtf"
CZ_GTF="${GENOME_DIR}/CZ_annotation.gtf"

# Prefixed pathogen FASTA and GTF files
CM_PREFIXED_FA="${GENOME_DIR}/CM_prefixed.fa"
CZ_PREFIXED_FA="${GENOME_DIR}/CZ_prefixed.fa"

CM_PREFIXED_GTF="${GENOME_DIR}/CM_prefixed.annotation.gtf"
CZ_PREFIXED_GTF="${GENOME_DIR}/CZ_prefixed.annotation.gtf"

# Combined host-pathogen references
ZM_CM_FA="${GENOME_DIR}/ZM_CM.fa"
ZM_CZ_FA="${GENOME_DIR}/ZM_CZ.fa"

# HISAT2 index prefixes
INDEX_ZM="${INDEX_DIR}/ZM"
INDEX_CM="${INDEX_DIR}/CM_prefixed"
INDEX_CZ="${INDEX_DIR}/CZ_prefixed"
INDEX_ZM_CM="${INDEX_DIR}/ZM_CM"
INDEX_ZM_CZ="${INDEX_DIR}/ZM_CZ"

# Trimmomatic adapter file
ADAPTER_FA="/path/to/trimmomatic/adapters/TruSeq3-PE-2.fa"

# Reference sequence lists for splitting BAM files
CM_REF_LIST="${GENOME_DIR}/cm_joint_refs.txt"
CZ_REF_LIST="${GENOME_DIR}/cz_joint_refs.txt"
MAIZE_FROM_CM_REF_LIST="${GENOME_DIR}/maize_from_CM_joint_refs.txt"
MAIZE_FROM_CZ_REF_LIST="${GENOME_DIR}/maize_from_CZ_joint_refs.txt"

# Output count files
MAIZE_COUNTS="${QUANT_DIR}/maize_counts.txt"
CM_COUNTS="${QUANT_DIR}/cm_counts.txt"
CZ_COUNTS="${QUANT_DIR}/cz_counts.txt"

# MultiQC report names
PRETRIM_MULTIQC_NAME="pretrim_report"
POSTTRIM_MULTIQC_NAME="posttrim_report"
