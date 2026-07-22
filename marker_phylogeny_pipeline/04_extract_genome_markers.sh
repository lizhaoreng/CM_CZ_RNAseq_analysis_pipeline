#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Extract ITS, TEF1, ACT, CAL, and HIS3 marker regions from CM/CZ genomes.
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
    [[ -s "${CM_GENOME}" ]] || die "Missing CM genome: ${CM_GENOME}"
    [[ -s "${CZ_GENOME}" ]] || die "Missing CZ genome: ${CZ_GENOME}"

    for marker in ${MARKERS_ALL}; do
        ref="${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta"
        [[ -s "${ref}" ]] || die "Missing marker reference: ${ref}"
    done
}

create_dirs() {
    mkdir -p genome_marker_blastdb
    mkdir -p genome_marker_blast
    mkdir -p genome_marker_summary
    mkdir -p "${GENOME_MARKER_DIR}/raw"
    mkdir -p "${GENOME_MARKER_DIR}/oriented"
    mkdir -p "${GENOME_MARKER_DIR}/renamed"
    mkdir -p genome_marker_validate
    mkdir -p marker_ref_blastdb
}

build_genome_blastdb() {
    log "Building genome BLAST databases..."

    makeblastdb -in "${CM_GENOME}" -dbtype nucl -out genome_marker_blastdb/CM_genome
    makeblastdb -in "${CZ_GENOME}" -dbtype nucl -out genome_marker_blastdb/CZ_genome
}

blast_markers_against_genomes() {
    log "BLAST marker references against genomes..."

    for marker in ${MARKERS_ALL}; do
        blastn \
            -query "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            -db genome_marker_blastdb/CM_genome \
            -out "genome_marker_blast/CM_${marker}_genome_blastn.tsv" \
            -evalue 1e-10 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -max_target_seqs 50 \
            -num_threads "${THREADS}"

        blastn \
            -query "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            -db genome_marker_blastdb/CZ_genome \
            -out "genome_marker_blast/CZ_${marker}_genome_blastn.tsv" \
            -evalue 1e-10 \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -max_target_seqs 50 \
            -num_threads "${THREADS}"
    done
}

select_best_hits() {
    log "Selecting best BLAST hit for each marker..."

    for sample in CM CZ; do
        for marker in ${MARKERS_ALL}; do
            infile="genome_marker_blast/${sample}_${marker}_genome_blastn.tsv"
            outfile="genome_marker_summary/${sample}_${marker}_besthit.tsv"

            if [[ -s "${infile}" ]]; then
                sort -k12,12nr "${infile}" | head -1 > "${outfile}"
            else
                : > "${outfile}"
                echo "WARNING: no BLAST hit for ${sample} ${marker}"
            fi
        done
    done
}

extract_regions() {
    log "Extracting genome-derived marker regions..."

    LOG="${GENOME_MARKER_DIR}/extraction_log.tsv"
    echo -e "sample\tmarker\tscaffold\tstart\tend\tstrand\tquery\tpident\tmatch_len\tbitscore" > "${LOG}"

    for sample in CM CZ; do
        if [[ "${sample}" == "CM" ]]; then
            genome="${CM_GENOME}"
        else
            genome="${CZ_GENOME}"
        fi

        for marker in ${MARKERS_ALL}; do
            blastfile="genome_marker_summary/${sample}_${marker}_besthit.tsv"

            if [[ ! -s "${blastfile}" ]]; then
                echo "WARNING: no best hit for ${sample} ${marker}, skip extraction"
                continue
            fi

            read -r qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore < "${blastfile}"

            if [[ "${sstart}" -le "${send}" ]]; then
                start="${sstart}"
                end="${send}"
                strand="+"
            else
                start="${send}"
                end="${sstart}"
                strand="-"
            fi

            start=$((start - FLANK))
            end=$((end + FLANK))

            if [[ "${start}" -lt 1 ]]; then
                start=1
            fi

            echo -e "${sample}\t${marker}\t${sseqid}\t${start}\t${end}\t${strand}\t${qseqid}\t${pident}\t${length}\t${bitscore}" \
                | tee -a "${LOG}"

            seqkit grep -p "${sseqid}" "${genome}" \
                | seqkit subseq -r "${start}:${end}" \
                > "${GENOME_MARKER_DIR}/raw/${sample}_genome_${marker}.raw.fasta"

            if [[ "${strand}" == "-" ]]; then
                seqkit seq -t DNA -r -p "${GENOME_MARKER_DIR}/raw/${sample}_genome_${marker}.raw.fasta" \
                    > "${GENOME_MARKER_DIR}/oriented/${sample}_genome_${marker}.fasta"
            else
                cp "${GENOME_MARKER_DIR}/raw/${sample}_genome_${marker}.raw.fasta" \
                   "${GENOME_MARKER_DIR}/oriented/${sample}_genome_${marker}.fasta"
            fi

            seqkit replace -p "^.*" -r "${sample}_genome_${marker}" \
                "${GENOME_MARKER_DIR}/oriented/${sample}_genome_${marker}.fasta" \
                > "${GENOME_MARKER_DIR}/renamed/${sample}_genome_${marker}.fasta"
        done
    done

    seqkit stats "${GENOME_MARKER_DIR}/renamed"/*.fasta
    grep "^>" "${GENOME_MARKER_DIR}/renamed"/*.fasta || true
}

validate_genome_markers() {
    log "Validating genome-derived markers against references..."

    for marker in ${MARKERS_ALL}; do
        makeblastdb \
            -in "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            -dbtype nucl \
            -out "marker_ref_blastdb/${marker}_ref"

        for sample in CM CZ; do
            query="${GENOME_MARKER_DIR}/renamed/${sample}_genome_${marker}.fasta"

            [[ -s "${query}" ]] || {
                echo "WARNING: missing query file: ${query}"
                continue
            }

            blastn \
                -query "${query}" \
                -db "marker_ref_blastdb/${marker}_ref" \
                -out "genome_marker_validate/${sample}_genome_${marker}_vs_ref.tsv" \
                -evalue 1e-10 \
                -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
        done
    done

    for marker in ${MARKERS_ALL}; do
        echo "===== CM genome ${marker} validation ====="
        sort -k12,12nr "genome_marker_validate/CM_genome_${marker}_vs_ref.tsv" | head -5 || true

        echo "===== CZ genome ${marker} validation ====="
        sort -k12,12nr "genome_marker_validate/CZ_genome_${marker}_vs_ref.tsv" | head -5 || true
    done
}

main() {
    check_dependencies
    check_inputs
    create_dirs
    build_genome_blastdb
    blast_markers_against_genomes
    select_best_hits
    extract_regions
    validate_genome_markers

    log "Genome marker extraction completed."
}

main "$@"
