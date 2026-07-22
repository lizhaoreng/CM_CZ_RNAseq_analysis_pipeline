#!/usr/bin/env bash

# ==============================================================================
# Configuration for multilocus marker extraction and phylogeny pipeline
# ==============================================================================

# Threads
THREADS=16
FASTQC_THREADS=32
FASTP_THREADS=16
TRINITY_THREADS=24
IQTREE_THREADS=AUTO

# Memory for Trinity
TRINITY_MAX_MEMORY="80G"

# Input in vitro RNA-seq reads
CM_R1_READS=("CM_1_R1.fastq.gz" "CM_2_R1.fastq.gz" "CM_3_R1.fastq.gz")
CM_R2_READS=("CM_1_R2.fastq.gz" "CM_2_R2.fastq.gz" "CM_3_R2.fastq.gz")

CZ_R1_READS=("CZ_1_R1.fastq.gz" "CZ_2_R1.fastq.gz" "CZ_3_R1.fastq.gz")
CZ_R2_READS=("CZ_1_R2.fastq.gz" "CZ_2_R2.fastq.gz" "CZ_3_R2.fastq.gz")

# Merged in vitro reads
CM_MERGED_R1="CM_invitro_all_R1.fastq.gz"
CM_MERGED_R2="CM_invitro_all_R2.fastq.gz"
CZ_MERGED_R1="CZ_invitro_all_R1.fastq.gz"
CZ_MERGED_R2="CZ_invitro_all_R2.fastq.gz"

# Clean reads
CM_CLEAN_R1="clean_reads/CM_R1.clean.fq.gz"
CM_CLEAN_R2="clean_reads/CM_R2.clean.fq.gz"
CZ_CLEAN_R1="clean_reads/CZ_R1.clean.fq.gz"
CZ_CLEAN_R2="clean_reads/CZ_R2.clean.fq.gz"

# Trinity output
CM_TRINITY="trinity_out/CM_trinity/Trinity.fasta"
CZ_TRINITY="trinity_out/CZ_trinity/Trinity.fasta"

# Reference genomes
CM_GENOME="CMCZgtf/CM_raw.fa"
CZ_GENOME="CMCZgtf/CZ_raw.fa"

# Markers
MARKERS_ALL="ITS TEF1 ACT CAL HIS3"
MARKERS_CODING="TEF1 ACT CAL HIS3"

# Flanking length for genome-derived marker extraction
FLANK=50

# Output directories
QC_RAW_DIR="qc_raw"
QC_CLEAN_DIR="qc_clean"
CLEAN_READS_DIR="clean_reads"
FASTP_REPORT_DIR="fastp_reports"
TRINITY_OUT_DIR="trinity_out"

# Entrez species list
SPECIES_LIST="species.list"

# Reference marker directories
SELECTED_DOWNLOADS_DIR="selected_downloads"
FINAL_MARKER_REF_DIR="final_marker_refs"

# Marker extraction directories
OWN_MARKER_DIR="own_marker_regions"
GENOME_MARKER_DIR="genome_marker_regions"

# Phylogeny output directories
FINAL_MARKER_WITH_OWN_DIR="final_marker_with_own"
ALIGNMENT_DIR="alignments"
TRIMMED_ALIGNMENT_DIR="alignments_trimmed"
CONCAT_DIR="concatenated"
IQTREE_DIR="iqtree_out"
CONCAT_NOCAL_DIR="concatenated_noCAL"
IQTREE_NOCAL_DIR="iqtree_noCAL"
