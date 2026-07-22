#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run InterProScan, KOfamScan, and eggNOG-mapper for fungal protein annotation
#
# Usage:
#   bash 01_run_interpro_kofam_eggnog.sh -g both -s all -t 64
#
# Options:
#   -g CM|CZ|both
#   -s interpro|kofam|eggnog|all
#   -t total threads
#   -f force overwrite
#   --no-parallel-genomes
# ==============================================================================

source ./config.sh

GENOME="both"
STEP="all"
FORCE=0
PARALLEL_GENOMES=1
LOCK_FILE="${ANNOTATION_DIR}/annotation.lock"

usage() {
    cat << EOF
Usage:
  bash 01_run_interpro_kofam_eggnog.sh [options]

Options:
  -g CM|CZ|both              Genome to process. Default: both
  -s interpro|kofam|eggnog|all
                              Annotation step. Default: all
  -t THREADS                 Total threads. Default from config.sh
  -f                         Force overwrite existing results
  --no-parallel-genomes      Do not run CM and CZ in parallel
  -h, --help                 Show help

Examples:
  bash 01_run_interpro_kofam_eggnog.sh -g CM -s interpro -t 32 -f
  bash 01_run_interpro_kofam_eggnog.sh -g both -s all -t 64
EOF
}

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -g) GENOME="$2"; shift 2 ;;
        -s) STEP="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -f) FORCE=1; shift ;;
        --no-parallel-genomes) PARALLEL_GENOMES=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown parameter: $1" ;;
    esac
done

exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
    die "Another annotation job is already running. Lock file: ${LOCK_FILE}"
fi

check_env() {
    log "Checking environment..."

    [[ -f "${CM_PROTEINS}" ]] || die "Missing CM protein file: ${CM_PROTEINS}"
    [[ -f "${CZ_PROTEINS}" ]] || die "Missing CZ protein file: ${CZ_PROTEINS}"
    [[ -x "${INTERPROSCAN}" ]] || die "InterProScan is not executable: ${INTERPROSCAN}"
    [[ -d "${KOFAM_DB}" ]] || die "KOFAM_DB not found: ${KOFAM_DB}"
    [[ -f "${KOFAM_DB}/ko_list" ]] || die "Missing KOfam ko_list"
    [[ -d "${KOFAM_DB}/profiles" ]] || die "Missing KOfam profiles"
    [[ -d "${EGGNOG_DATA_DIR}" ]] || die "EGGNOG_DATA_DIR not found"

    command -v emapper.py >/dev/null 2>&1 || die "emapper.py not found"
    command -v diamond >/dev/null 2>&1 || die "diamond not found"

    log "InterProScan version:"
    "${INTERPROSCAN}" -version || true

    log "Java version:"
    java -version 2>&1 | head -n 2 || true

    log "Total threads: ${THREADS}"
}

protein_file_for_genome() {
    local g="$1"
    if [[ "${g}" == "CM" ]]; then
        echo "${CM_PROTEINS}"
    elif [[ "${g}" == "CZ" ]]; then
        echo "${CZ_PROTEINS}"
    else
        die "Unknown genome: ${g}"
    fi
}

thread_plan() {
    local genome_count="$1"
    if [[ "${genome_count}" -eq 2 && "${PARALLEL_GENOMES}" -eq 1 ]]; then
        echo $(( THREADS / 2 ))
    else
        echo "${THREADS}"
    fi
}

run_cmd() {
    local name="$1"
    local cmd="$2"
    local log_file="$3"

    log "Start: ${name}"
    log "Log file: ${log_file}"
    echo "${cmd}" > "${log_file}.cmd"

    bash -c "${cmd}" > "${log_file}" 2>&1

    log "Finished: ${name}"
}

run_interpro() {
    local g="$1"
    local cpu="$2"
    local input
    input=$(protein_file_for_genome "${g}")

    local prefix="${INTERPRO_DIR}/${g}_interpro"
    local tmp="${TMP_BASE}/interpro_${g}_$$"
    local log_file="${LOG_DIR}/${g}.interpro.log"

    if [[ -s "${prefix}.tsv" && "${FORCE}" -eq 0 ]]; then
        log "${g} InterProScan result exists. Skipping: ${prefix}.tsv"
        return 0
    fi

    rm -rf "${tmp}"
    mkdir -p "${tmp}"

    local cmd="${INTERPROSCAN} \
        -i \"${input}\" \
        -f TSV,GFF3 \
        -appl Pfam,TIGRFAM,ProSitePatterns \
        -d \"${INTERPRO_DIR}\" \
        -T \"${tmp}\" \
        -cpu \"${cpu}\" \
        -goterms \
        -iprlookup"

    run_cmd "${g} InterProScan" "${cmd}" "${log_file}"

    local base
    base=$(basename "${input}")

    [[ -f "${INTERPRO_DIR}/${base}.tsv" ]] && mv -f "${INTERPRO_DIR}/${base}.tsv" "${prefix}.tsv"
    [[ -f "${INTERPRO_DIR}/${base}.gff3" ]] && mv -f "${INTERPRO_DIR}/${base}.gff3" "${prefix}.gff3"

    rm -rf "${tmp}"

    [[ -s "${prefix}.tsv" ]] || die "${g} InterProScan result is missing or empty"
}

run_kofam() {
    local g="$1"
    local cpu="$2"
    local input
    input=$(protein_file_for_genome "${g}")

    local output="${KOFAM_DIR}/${g}_kofam.tsv"
    local tmp="${TMP_BASE}/kofam_${g}_$$"
    local log_file="${LOG_DIR}/${g}.kofam.log"

    if [[ -s "${output}" && "${FORCE}" -eq 0 ]]; then
        log "${g} KOfamScan result exists. Skipping: ${output}"
        return 0
    fi

    rm -rf "${tmp}"
    mkdir -p "${tmp}"

    local cmd="exec_annotation \
        --profile=\"${KOFAM_DB}/profiles\" \
        --ko-list=\"${KOFAM_DB}/ko_list\" \
        -f detail-tsv \
        --cpu \"${cpu}\" \
        --tmp-dir=\"${tmp}\" \
        -o \"${output}\" \
        \"${input}\""

    run_cmd "${g} KOfamScan" "${cmd}" "${log_file}"

    rm -rf "${tmp}"

    [[ -s "${output}" ]] || die "${g} KOfamScan result is missing or empty"
}

run_eggnog() {
    local g="$1"
    local cpu="$2"
    local input
    input=$(protein_file_for_genome "${g}")

    local out_prefix="${g}_eggnog"
    local final="${EGGNOG_DIR}/${g}_eggnog.emapper.annotations"
    local tmp="${TMP_BASE}/eggnog_${g}_$$"
    local temp_out="${EGGNOG_DIR}/${g}_work"
    local log_file="${LOG_DIR}/${g}.eggnog.log"

    if [[ -s "${final}" && "${FORCE}" -eq 0 ]]; then
        log "${g} eggNOG result exists. Skipping: ${final}"
        return 0
    fi

    rm -rf "${tmp}" "${temp_out}"
    mkdir -p "${tmp}" "${temp_out}"

    local pfam_opt=""
    if [[ -d "${EGGNOG_DATA_DIR}/pfam" ]]; then
        pfam_opt="--pfam_realign realign"
    fi

    local cmd="emapper.py \
        -i \"${input}\" \
        -o \"${out_prefix}\" \
        -m diamond \
        --cpu \"${cpu}\" \
        --tax_scope 4751 \
        --go_evidence non-electronic \
        --target_orthologs all \
        --seed_ortholog_evalue 1e-5 \
        --seed_ortholog_score 60 \
        --sensmode fast \
        --block_size 4.0 \
        --pident 50 \
        --query_cover 60 \
        --subject_cover 60 \
        --temp_dir \"${tmp}\" \
        --report_orthologs \
        --no_file_comments \
        --excel \
        --override \
        --output_dir \"${temp_out}\" \
        --data_dir \"${EGGNOG_DATA_DIR}\" \
        ${pfam_opt}"

    run_cmd "${g} eggNOG-mapper" "${cmd}" "${log_file}"

    mv -f "${temp_out}/${out_prefix}.emapper.annotations" "${EGGNOG_DIR}/" 2>/dev/null || true
    mv -f "${temp_out}/${out_prefix}.emapper.seed_orthologs" "${EGGNOG_DIR}/" 2>/dev/null || true
    mv -f "${temp_out}/${out_prefix}.emapper.hits" "${EGGNOG_DIR}/" 2>/dev/null || true
    mv -f "${temp_out}/${out_prefix}.emapper.xlsx" "${EGGNOG_DIR}/" 2>/dev/null || true

    rm -rf "${tmp}" "${temp_out}"

    [[ -s "${final}" ]] || die "${g} eggNOG result is missing or empty"
}

process_one_genome() {
    local g="$1"
    local cpu="$2"
    local input
    input=$(protein_file_for_genome "${g}")

    local n
    n=$(grep -c '^>' "${input}" || true)

    log "=============================="
    log "Genome: ${g}"
    log "Protein number: ${n}"
    log "Threads: ${cpu}"
    log "Step: ${STEP}"
    log "=============================="

    case "${STEP}" in
        interpro) run_interpro "${g}" "${cpu}" ;;
        kofam) run_kofam "${g}" "${cpu}" ;;
        eggnog) run_eggnog "${g}" "${cpu}" ;;
        all)
            run_interpro "${g}" "${cpu}"
            run_kofam "${g}" "${cpu}"
            run_eggnog "${g}" "${cpu}"
            ;;
        *) die "Unknown step: ${STEP}" ;;
    esac

    log "${g} completed"
}

summary() {
    log "Result summary:"

    for g in CM CZ; do
        [[ "${GENOME}" == "${g}" || "${GENOME}" == "both" ]] || continue

        echo "---- ${g} ----"
        for f in \
            "${INTERPRO_DIR}/${g}_interpro.tsv" \
            "${KOFAM_DIR}/${g}_kofam.tsv" \
            "${EGGNOG_DIR}/${g}_eggnog.emapper.annotations"
        do
            if [[ -s "${f}" ]]; then
                echo "$(wc -l < "${f}") lines  ${f}"
            else
                echo "missing/empty  ${f}"
            fi
        done
    done
}

main() {
    check_env

    local genomes=()
    case "${GENOME}" in
        CM) genomes=("CM") ;;
        CZ) genomes=("CZ") ;;
        both) genomes=("CM" "CZ") ;;
        *) die "Unknown genome: ${GENOME}" ;;
    esac

    local genome_count="${#genomes[@]}"
    local cpu_per_genome
    cpu_per_genome=$(thread_plan "${genome_count}")

    if [[ "${cpu_per_genome}" -lt 1 ]]; then
        cpu_per_genome=1
    fi

    log "Genome number: ${genome_count}"
    log "Threads per genome: ${cpu_per_genome}"
    log "Parallel genomes: ${PARALLEL_GENOMES}"

    if [[ "${genome_count}" -eq 2 && "${PARALLEL_GENOMES}" -eq 1 ]]; then
        process_one_genome "CM" "${cpu_per_genome}" > "${LOG_DIR}/CM.main.log" 2>&1 &
        pid_cm=$!

        process_one_genome "CZ" "${cpu_per_genome}" > "${LOG_DIR}/CZ.main.log" 2>&1 &
        pid_cz=$!

        wait "${pid_cm}"
        wait "${pid_cz}"
    else
        for g in "${genomes[@]}"; do
            process_one_genome "${g}" "${cpu_per_genome}"
        done
    fi

    summary
    log "All annotation jobs completed"
}

main "$@"
