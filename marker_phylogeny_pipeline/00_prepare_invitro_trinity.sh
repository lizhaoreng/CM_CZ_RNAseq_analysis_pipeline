#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Prepare in vitro RNA-seq reads and run Trinity de novo transcriptome assembly
#
# Steps:
#   1. Merge biological replicate reads for each species
#   2. Run FastQC/MultiQC before trimming
#   3. Trim reads using fastp
#   4. Run FastQC/MultiQC after trimming
#   5. Run Trinity assembly for CM and CZ
#   6. Generate Trinity assembly statistics
#   7. Extract longest isoform per Trinity gene
#   8. Build BLAST nucleotide databases
# ==============================================================================

source ./config.sh

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || die "$1 not found"
}

check_inputs() {
    log "Checking required tools..."

    for cmd in fastqc multiqc fastp Trinity makeblastdb; do
        check_command "$cmd"
    done

    if [[ -z "${TRINITY_HOME:-}" ]]; then
        die "TRINITY_HOME is not set. Please export TRINITY_HOME before running this script."
    fi

    [[ -x "${TRINITY_HOME}/util/TrinityStats.pl" ]] || die "TrinityStats.pl not found"
    [[ -x "${TRINITY_HOME}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl" ]] || die "get_longest_isoform script not found"

    log "Checking input FASTQ files..."

    for f in "${CM_R1_READS[@]}" "${CM_R2_READS[@]}" "${CZ_R1_READS[@]}" "${CZ_R2_READS[@]}"; do
        [[ -s "$f" ]] || die "Missing FASTQ file: $f"
    done
}

merge_reads() {
    log "Merging in vitro RNA-seq reads..."

    cat "${CM_R1_READS[@]}" > "${CM_MERGED_R1}"
    cat "${CM_R2_READS[@]}" > "${CM_MERGED_R2}"

    cat "${CZ_R1_READS[@]}" > "${CZ_MERGED_R1}"
    cat "${CZ_R2_READS[@]}" > "${CZ_MERGED_R2}"

    log "Merged read files generated:"
    ls -lh "${CM_MERGED_R1}" "${CM_MERGED_R2}" "${CZ_MERGED_R1}" "${CZ_MERGED_R2}"
}

run_fastqc_raw() {
    log "Running FastQC on raw merged reads..."

    mkdir -p "${QC_RAW_DIR}"

    fastqc \
        "${CM_MERGED_R1}" "${CM_MERGED_R2}" \
        "${CZ_MERGED_R1}" "${CZ_MERGED_R2}" \
        -o "${QC_RAW_DIR}" \
        -t "${FASTQC_THREADS}"

    multiqc "${QC_RAW_DIR}" -o "${QC_RAW_DIR}" --force
}

run_fastp() {
    log "Running fastp trimming..."

    mkdir -p "${CLEAN_READS_DIR}" "${FASTP_REPORT_DIR}"

    fastp \
        -i "${CM_MERGED_R1}" \
        -I "${CM_MERGED_R2}" \
        -o "${CM_CLEAN_R1}" \
        -O "${CM_CLEAN_R2}" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality 20 \
        --length_required 50 \
        --thread "${FASTP_THREADS}" \
        --html "${FASTP_REPORT_DIR}/CM_fastp.html" \
        --json "${FASTP_REPORT_DIR}/CM_fastp.json"

    fastp \
        -i "${CZ_MERGED_R1}" \
        -I "${CZ_MERGED_R2}" \
        -o "${CZ_CLEAN_R1}" \
        -O "${CZ_CLEAN_R2}" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality 20 \
        --length_required 50 \
        --thread "${FASTP_THREADS}" \
        --html "${FASTP_REPORT_DIR}/CZ_fastp.html" \
        --json "${FASTP_REPORT_DIR}/CZ_fastp.json"
}

run_fastqc_clean() {
    log "Running FastQC on cleaned reads..."

    mkdir -p "${QC_CLEAN_DIR}"

    fastqc \
        "${CM_CLEAN_R1}" "${CM_CLEAN_R2}" \
        "${CZ_CLEAN_R1}" "${CZ_CLEAN_R2}" \
        -o "${QC_CLEAN_DIR}" \
        -t "${FASTQC_THREADS}"

    multiqc "${QC_CLEAN_DIR}" "${FASTP_REPORT_DIR}" -o "${QC_CLEAN_DIR}" --force
}

run_trinity() {
    log "Running Trinity assembly..."

    mkdir -p "${TRINITY_OUT_DIR}"

    Trinity \
        --seqType fq \
        --left "${CM_CLEAN_R1}" \
        --right "${CM_CLEAN_R2}" \
        --CPU "${TRINITY_THREADS}" \
        --max_memory "${TRINITY_MAX_MEMORY}" \
        --output "${TRINITY_OUT_DIR}/CM_trinity"

    Trinity \
        --seqType fq \
        --left "${CZ_CLEAN_R1}" \
        --right "${CZ_CLEAN_R2}" \
        --CPU "${TRINITY_THREADS}" \
        --max_memory "${TRINITY_MAX_MEMORY}" \
        --output "${TRINITY_OUT_DIR}/CZ_trinity"
}

trinity_stats_and_longest_isoform() {
    log "Generating Trinity statistics..."

    "${TRINITY_HOME}/util/TrinityStats.pl" "${CM_TRINITY}" > CM_trinity_stats.txt
    "${TRINITY_HOME}/util/TrinityStats.pl" "${CZ_TRINITY}" > CZ_trinity_stats.txt

    log "Extracting longest isoform per Trinity gene..."

    "${TRINITY_HOME}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl" \
        "${CM_TRINITY}" > "${TRINITY_OUT_DIR}/CM_trinity.longest.fasta"

    "${TRINITY_HOME}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl" \
        "${CZ_TRINITY}" > "${TRINITY_OUT_DIR}/CZ_trinity.longest.fasta"
}

build_trinity_blastdb() {
    log "Building BLAST databases for Trinity assemblies..."

    mkdir -p "${TRINITY_OUT_DIR}/blastdb"

    makeblastdb \
        -in "${TRINITY_OUT_DIR}/CM_trinity.longest.fasta" \
        -dbtype nucl \
        -out "${TRINITY_OUT_DIR}/blastdb/CM_trinity_longest"

    makeblastdb \
        -in "${TRINITY_OUT_DIR}/CZ_trinity.longest.fasta" \
        -dbtype nucl \
        -out "${TRINITY_OUT_DIR}/blastdb/CZ_trinity_longest"
}

main() {
    check_inputs
    merge_reads
    run_fastqc_raw
    run_fastp
    run_fastqc_clean
    run_trinity
    trinity_stats_and_longest_isoform
    build_trinity_blastdb

    log "In vitro Trinity preparation completed."
}

main "$@"
