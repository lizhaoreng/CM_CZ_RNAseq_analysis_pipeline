#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run DIAMOND BLASTP against NCBI RefSeq fungi protein database
# ==============================================================================

source ./config.sh

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_inputs() {
    command -v diamond >/dev/null 2>&1 || die "diamond not found"

    [[ -f "${NCBI_FUNGI_DMND}" ]] || die "NCBI fungi DIAMOND database not found: ${NCBI_FUNGI_DMND}"
    [[ -f "${CM_PROTEINS}" ]] || die "CM protein file not found: ${CM_PROTEINS}"
    [[ -f "${CZ_PROTEINS}" ]] || die "CZ protein file not found: ${CZ_PROTEINS}"

    mkdir -p "${BLAST_DIR}"
}

run_blast_one() {
    local species="$1"
    local query="$2"

    local base
    base=$(basename "${query}" .faa)

    local raw_out="${BLAST_DIR}/${base}_blast_results.txt"
    local header_out="${BLAST_DIR}/${base}_blast_with_header.txt"
    local high_sim_out="${BLAST_DIR}/${base}_high_similarity.txt"
    local log_file="${LOG_DIR}/${base}.ncbi_blast.log"

    log "Running NCBI fungi DIAMOND BLAST for ${species}"
    log "Query: ${query}"

    diamond blastp \
        --query "${query}" \
        --db "${NCBI_FUNGI_DMND}" \
        --out "${raw_out}" \
        --evalue "${EVALUE}" \
        --threads "${THREADS}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --max-target-seqs "${MAX_TARGET_SEQS_NCBI}" \
        --more-sensitive \
        --verbose \
        > "${log_file}" 2>&1

    echo -e "Query_ID\tSubject_ID\tIdentity%\tLength\tMismatches\tGaps\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tDescription" > "${header_out}"
    cat "${raw_out}" >> "${header_out}"

    echo -e "Query_ID\tSubject_ID\tIdentity%\tLength\tMismatches\tGaps\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tDescription" > "${high_sim_out}"
    awk -F'\t' '$3 > 70' "${raw_out}" >> "${high_sim_out}"

    local query_count
    local hit_count
    local unique_queries

    query_count=$(grep -c '^>' "${query}" || true)
    hit_count=$(wc -l < "${raw_out}")
    unique_queries=$(cut -f1 "${raw_out}" | sort -u | wc -l)

    log "${species} query proteins: ${query_count}"
    log "${species} total hits: ${hit_count}"
    log "${species} queries with hits: ${unique_queries}"
    log "Output: ${header_out}"
}

write_summary() {
    local summary="${BLAST_DIR}/blast_summary_report.txt"

    {
        echo "NCBI fungi DIAMOND BLAST summary"
        echo "================================"
        echo "Date: $(date)"
        echo "Database: ${NCBI_FUNGI_DMND}"
        echo "E-value: ${EVALUE}"
        echo "Max target sequences: ${MAX_TARGET_SEQS_NCBI}"
        echo "Threads: ${THREADS}"
        echo ""

        for f in "${BLAST_DIR}"/*_blast_with_header.txt; do
            [[ -f "${f}" ]] || continue
            local base
            base=$(basename "${f}" _blast_with_header.txt)
            local hits
            local unique
            hits=$(tail -n +2 "${f}" | wc -l)
            unique=$(tail -n +2 "${f}" | cut -f1 | sort -u | wc -l)

            echo "${base}"
            echo "  Total hits: ${hits}"
            echo "  Queries with hits: ${unique}"
            echo ""
        done
    } > "${summary}"

    log "Summary report saved: ${summary}"
}

main() {
    check_inputs

    run_blast_one "CM" "${CM_PROTEINS}"
    run_blast_one "CZ" "${CZ_PROTEINS}"

    write_summary

    log "NCBI fungi DIAMOND BLAST completed"
}

main "$@"
