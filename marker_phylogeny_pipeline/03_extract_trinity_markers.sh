#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Extract TEF1, ACT, CAL, and HIS3 marker regions from CM/CZ Trinity assemblies.
#
# Note:
#   The final extraction coordinates below are manually curated based on BLAST
#   results against reference marker sequences.
# ==============================================================================

source ./config.sh

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_dependencies() {
    for cmd in makeblastdb blastn seqkit; do
        command -v "$cmd" >/dev/null 2>&1 || die "$cmd not found"
    done
}

check_inputs() {
    [[ -s "${CM_TRINITY}" ]] || die "Missing CM Trinity assembly: ${CM_TRINITY}"
    [[ -s "${CZ_TRINITY}" ]] || die "Missing CZ Trinity assembly: ${CZ_TRINITY}"

    for marker in ${MARKERS_CODING}; do
        ref="${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta"
        [[ -s "${ref}" ]] || die "Missing marker reference: ${ref}"
    done
}

build_blastdb() {
    log "Building BLAST databases for Trinity assemblies..."

    mkdir -p blastdb marker_blast

    makeblastdb -in "${CM_TRINITY}" -dbtype nucl -out blastdb/CM_trinity
    makeblastdb -in "${CZ_TRINITY}" -dbtype nucl -out blastdb/CZ_trinity
}

blast_markers_against_trinity() {
    log "BLAST marker references against Trinity assemblies..."

    for marker in ${MARKERS_CODING}; do
        blastn \
            -query "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            -db blastdb/CM_trinity \
            -out "marker_blast/CM_${marker}_blastn.tsv" \
            -evalue 1e-20 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -max_target_seqs 20 \
            -num_threads "${THREADS}"

        blastn \
            -query "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            -db blastdb/CZ_trinity \
            -out "marker_blast/CZ_${marker}_blastn.tsv" \
            -evalue 1e-20 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -max_target_seqs 20 \
            -num_threads "${THREADS}"
    done
}

show_top_hits() {
    log "Showing top Trinity BLAST hits..."

    for marker in ${MARKERS_CODING}; do
        echo "===== CM ${marker} ====="
        sort -k12,12nr "marker_blast/CM_${marker}_blastn.tsv" | head -10 || true

        echo "===== CZ ${marker} ====="
        sort -k12,12nr "marker_blast/CZ_${marker}_blastn.tsv" | head -10 || true
    done
}

extract_confirmed_regions() {
    log "Extracting manually curated marker regions..."

    mkdir -p "${OWN_MARKER_DIR}/raw" "${OWN_MARKER_DIR}/oriented" "${OWN_MARKER_DIR}/renamed"

    # --------------------------------------------------------------------------
    # TEF1
    # --------------------------------------------------------------------------
    seqkit grep -p "TRINITY_DN483_c1_g1_i118" "${CM_TRINITY}" \
        | seqkit subseq -r 298:608 \
        > "${OWN_MARKER_DIR}/raw/CM_TEF1.raw.fasta"

    seqkit grep -p "TRINITY_DN16329_c0_g1_i1" "${CZ_TRINITY}" \
        | seqkit subseq -r 27:334 \
        > "${OWN_MARKER_DIR}/raw/CZ_TEF1.raw.fasta"

    # --------------------------------------------------------------------------
    # ACT
    # --------------------------------------------------------------------------
    seqkit grep -p "TRINITY_DN699_c0_g1_i2" "${CM_TRINITY}" \
        | seqkit subseq -r 1291:1516 \
        > "${OWN_MARKER_DIR}/raw/CM_ACT.raw.fasta"

    seqkit grep -p "TRINITY_DN1613_c2_g1_i49" "${CZ_TRINITY}" \
        | seqkit subseq -r 3904:4126 \
        > "${OWN_MARKER_DIR}/raw/CZ_ACT.raw.fasta"

    # --------------------------------------------------------------------------
    # CAL
    # --------------------------------------------------------------------------
    seqkit grep -p "TRINITY_DN2150_c0_g1_i1" "${CM_TRINITY}" \
        | seqkit subseq -r 4295:4603 \
        > "${OWN_MARKER_DIR}/raw/CM_CAL.raw.fasta"

    seqkit grep -p "TRINITY_DN3027_c0_g1_i1" "${CZ_TRINITY}" \
        | seqkit subseq -r 3440:3748 \
        > "${OWN_MARKER_DIR}/raw/CZ_CAL.raw.fasta"

    seqkit grep -p "TRINITY_DN3682_c0_g1_i2" "${CM_TRINITY}" \
        | seqkit subseq -r 1:141 \
        > "${OWN_MARKER_DIR}/raw/CM_CAL_short.raw.fasta"

    # --------------------------------------------------------------------------
    # HIS3
    # --------------------------------------------------------------------------
    seqkit grep -p "TRINITY_DN274_c6_g1_i1" "${CM_TRINITY}" \
        | seqkit subseq -r 4515:4747 \
        > "${OWN_MARKER_DIR}/raw/CM_HIS3.raw.fasta"

    seqkit grep -p "TRINITY_DN3431_c0_g1_i2" "${CZ_TRINITY}" \
        | seqkit subseq -r 676:908 \
        > "${OWN_MARKER_DIR}/raw/CZ_HIS3.raw.fasta"

    seqkit stats "${OWN_MARKER_DIR}/raw"/*.fasta
}

orient_sequences() {
    log "Orienting marker sequences..."

    cp "${OWN_MARKER_DIR}/raw/CM_TEF1.raw.fasta" "${OWN_MARKER_DIR}/oriented/CM_TEF1.fasta"
    cp "${OWN_MARKER_DIR}/raw/CZ_ACT.raw.fasta"  "${OWN_MARKER_DIR}/oriented/CZ_ACT.fasta"
    cp "${OWN_MARKER_DIR}/raw/CM_CAL.raw.fasta"  "${OWN_MARKER_DIR}/oriented/CM_CAL.fasta"

    seqkit seq -t DNA -r -p "${OWN_MARKER_DIR}/raw/CZ_TEF1.raw.fasta" \
        > "${OWN_MARKER_DIR}/oriented/CZ_TEF1.fasta"

    seqkit seq -t DNA -r -p "${OWN_MARKER_DIR}/raw/CM_ACT.raw.fasta" \
        > "${OWN_MARKER_DIR}/oriented/CM_ACT.fasta"

    seqkit seq -t DNA -r -p "${OWN_MARKER_DIR}/raw/CZ_CAL.raw.fasta" \
        > "${OWN_MARKER_DIR}/oriented/CZ_CAL.fasta"

    seqkit seq -t DNA -r -p "${OWN_MARKER_DIR}/raw/CM_HIS3.raw.fasta" \
        > "${OWN_MARKER_DIR}/oriented/CM_HIS3.fasta"

    seqkit seq -t DNA -r -p "${OWN_MARKER_DIR}/raw/CZ_HIS3.raw.fasta" \
        > "${OWN_MARKER_DIR}/oriented/CZ_HIS3.fasta"

    seqkit seq -t DNA -r -p "${OWN_MARKER_DIR}/raw/CM_CAL_short.raw.fasta" \
        > "${OWN_MARKER_DIR}/oriented/CM_CAL_short.fasta"
}

rename_sequences() {
    log "Renaming Trinity-derived marker sequences..."

    for marker in ${MARKERS_CODING}; do
        seqkit replace -p "^.*" -r "CM_${marker}" \
            "${OWN_MARKER_DIR}/oriented/CM_${marker}.fasta" \
            > "${OWN_MARKER_DIR}/renamed/CM_${marker}.fasta"

        seqkit replace -p "^.*" -r "CZ_${marker}" \
            "${OWN_MARKER_DIR}/oriented/CZ_${marker}.fasta" \
            > "${OWN_MARKER_DIR}/renamed/CZ_${marker}.fasta"
    done

    seqkit replace -p "^.*" -r "CM_CAL_short" \
        "${OWN_MARKER_DIR}/oriented/CM_CAL_short.fasta" \
        > "${OWN_MARKER_DIR}/renamed/CM_CAL_short.fasta"

    seqkit stats "${OWN_MARKER_DIR}/renamed"/*.fasta
    grep "^>" "${OWN_MARKER_DIR}/renamed"/*.fasta || true
}

validate_trinity_markers() {
    log "Validating Trinity-derived markers against references..."

    mkdir -p marker_ref_blastdb marker_validate

    for marker in ${MARKERS_CODING}; do
        makeblastdb \
            -in "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            -dbtype nucl \
            -out "marker_ref_blastdb/${marker}_ref"

        blastn \
            -query "${OWN_MARKER_DIR}/renamed/CM_${marker}.fasta" \
            -db "marker_ref_blastdb/${marker}_ref" \
            -out "marker_validate/CM_${marker}_vs_ref.tsv" \
            -evalue 1e-10 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"

        blastn \
            -query "${OWN_MARKER_DIR}/renamed/CZ_${marker}.fasta" \
            -db "marker_ref_blastdb/${marker}_ref" \
            -out "marker_validate/CZ_${marker}_vs_ref.tsv" \
            -evalue 1e-10 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
    done

    for marker in ${MARKERS_CODING}; do
        echo "===== CM ${marker} validation ====="
        sort -k12,12nr "marker_validate/CM_${marker}_vs_ref.tsv" | head -5 || true

        echo "===== CZ ${marker} validation ====="
        sort -k12,12nr "marker_validate/CZ_${marker}_vs_ref.tsv" | head -5 || true
    done
}

main() {
    check_dependencies
    check_inputs
    build_blastdb
    blast_markers_against_trinity
    show_top_hits
    extract_confirmed_regions
    orient_sequences
    rename_sequences
    validate_trinity_markers

    log "Trinity marker extraction completed."
}

main "$@"
