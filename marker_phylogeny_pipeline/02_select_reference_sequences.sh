#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Select curated reference marker sequences from downloaded NCBI records.
#
# Input:
#   selected_downloads/{ITS,TEF1,ACT,CAL,HIS3}/*.fasta
#
# Output:
#   marker_tables/
#   screening_tables/
#   selected_refs_by_acc/
#   selected_marker_refs/
#   selected_other_refs/
#   final_marker_refs/
# ==============================================================================

source ./config.sh

MARKERS="TEF1 ACT CAL HIS3 ITS"

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_dependencies() {
    command -v seqkit >/dev/null 2>&1 || die "seqkit not found"
}

extract_headers() {
    log "Extracting FASTA headers..."

    mkdir -p headers marker_tables

    for marker in ${MARKERS}; do
        mkdir -p "headers/${marker}"
        echo -e "marker\tfile\taccession\theader" > "marker_tables/${marker}_headers.tsv"

        for f in "${SELECTED_DOWNLOADS_DIR}/${marker}"/*.fasta; do
            [[ -e "${f}" ]] || continue

            base=$(basename "${f}" .fasta)
            grep "^>" "${f}" > "headers/${marker}/${base}.headers.txt" || true

            grep "^>" "${f}" | sed 's/^>//' | awk -v m="${marker}" -v file="$(basename "${f}")" '{
                split($1,a,".");
                print m "\t" file "\t" a[1] "\t" $0
            }' >> "marker_tables/${marker}_headers.tsv" || true
        done
    done
}

build_candidate_tables() {
    log "Building candidate tables..."

    mkdir -p screening_tables screening_tables/filtered

    for marker in ${MARKERS}; do
        echo -e "marker\tfile\taccession\tlength\theader" > "screening_tables/${marker}_candidates.tsv"

        for f in "${SELECTED_DOWNLOADS_DIR}/${marker}"/*.fasta; do
            [[ -e "${f}" ]] || continue

            base=$(basename "${f}")

            seqkit fx2tab -n -l -i "${f}" | while read -r id len; do
                header=$(grep -m 1 "^>${id}" "${f}" | sed 's/^>//' || true)
                acc=$(echo "${id}" | sed 's/\..*//')
                echo -e "${marker}\t${base}\t${acc}\t${len}\t${header}" >> "screening_tables/${marker}_candidates.tsv"
            done
        done
    done
}

filter_candidates() {
    log "Filtering candidate sequences..."

    awk -F'\t' 'NR==1 || ($5 !~ /UNVERIFIED/ && $4 >= 450)' screening_tables/ITS_candidates.tsv  > screening_tables/filtered/ITS_filtered.tsv
    awk -F'\t' 'NR==1 || ($5 !~ /UNVERIFIED/ && $4 >= 250)' screening_tables/TEF1_candidates.tsv > screening_tables/filtered/TEF1_filtered.tsv
    awk -F'\t' 'NR==1 || ($5 !~ /UNVERIFIED/ && $4 >= 200)' screening_tables/ACT_candidates.tsv  > screening_tables/filtered/ACT_filtered.tsv
    awk -F'\t' 'NR==1 || ($5 !~ /UNVERIFIED/ && $4 >= 250)' screening_tables/CAL_candidates.tsv  > screening_tables/filtered/CAL_filtered.tsv
    awk -F'\t' 'NR==1 || ($5 !~ /UNVERIFIED/ && $4 >= 250)' screening_tables/HIS3_candidates.tsv > screening_tables/filtered/HIS3_filtered.tsv
}

prepare_core_accessions() {
    log "Preparing core accession list..."

    mkdir -p accession_lists selected_refs_by_acc all_marker_pool selected_marker_refs

    cat > accession_lists/core_refs.acc <<'ACC'
OR945717
PV133732
DQ185080
DQ185092
DQ185104
DQ185128
OQ773855
KM087697
KM087649
KM087673
NR_111205
DQ185093
DQ185105
DQ185117
DQ185129
EU569227
EU569218
DQ185094
DQ185106
DQ185118
DQ185130
ACC

    cat accession_lists/core_refs.acc
}

build_all_marker_pool() {
    log "Building all-marker sequence pool..."

    cat \
        "${SELECTED_DOWNLOADS_DIR}/TEF1"/*.fasta \
        "${SELECTED_DOWNLOADS_DIR}/ACT"/*.fasta \
        "${SELECTED_DOWNLOADS_DIR}/CAL"/*.fasta \
        "${SELECTED_DOWNLOADS_DIR}/HIS3"/*.fasta \
        "${SELECTED_DOWNLOADS_DIR}/ITS"/*.fasta \
        > all_marker_pool/all_downloaded_markers.fasta
}

extract_core_refs() {
    log "Extracting core references by accession..."

    seqkit grep -r -f accession_lists/core_refs.acc all_marker_pool/all_downloaded_markers.fasta \
        > selected_refs_by_acc/core_refs.fasta

    log "Core reference count:"
    grep -c "^>" selected_refs_by_acc/core_refs.fasta || true
}

split_core_refs_by_marker() {
    log "Splitting core references by marker..."

    seqkit grep -n -r -i -p "elongation factor|tef1|EF1" selected_refs_by_acc/core_refs.fasta \
        > selected_marker_refs/TEF1_ref.selected.fasta

    seqkit grep -n -r -i -p "actin|act\)" selected_refs_by_acc/core_refs.fasta \
        > selected_marker_refs/ACT_ref.selected.fasta

    seqkit grep -n -r -i -p "calmodulin|cmdA|CAL" selected_refs_by_acc/core_refs.fasta \
        > selected_marker_refs/CAL_ref.selected.fasta

    seqkit grep -n -r -i -p "histone H3|HIS" selected_refs_by_acc/core_refs.fasta \
        > selected_marker_refs/HIS3_ref.selected.fasta

    seqkit grep -n -r -i -p "internal transcribed spacer|ITS region|ribosomal RNA" selected_refs_by_acc/core_refs.fasta \
        > selected_marker_refs/ITS_ref.selected.fasta
}

select_other_species_refs() {
    log "Selecting additional species references..."

    cat > other_species.list <<'LIST'
Cercospora beticola
Cercospora kikuchii
Cercospora sojina
Cercospora apii
Cercospora asparagi
Cercospora flagellaris
Cercospora cf. flagellaris
Cercospora nicotianae
Pseudocercospora fijiensis
Pseudocercospora griseola
LIST

    mkdir -p selected_other_refs/{TEF1,ACT,CAL,HIS3,ITS}

    echo -e "species\tmarker\tn_selected\toutfile" > selected_other_refs/selection_log.tsv

    while read -r sp; do
        [[ -z "${sp}" ]] && continue

        sp_file=$(echo "${sp}" | sed 's/ /_/g; s/\./_/g')
        sp_pattern=$(echo "${sp}" | sed 's/\./\\./g')

        for marker in TEF1 ACT CAL HIS3 ITS; do
            infile="${SELECTED_DOWNLOADS_DIR}/${marker}/${sp_file}_${marker}.fasta"
            outfile="selected_other_refs/${marker}/${sp_file}_${marker}.fasta"

            if [[ -s "${infile}" ]]; then
                seqkit grep -n -r -i -p "${sp_pattern}" "${infile}" \
                    | seqkit grep -n -r -v -i -p "UNVERIFIED" \
                    | seqkit seq -m 200 \
                    | seqkit head -n 1 > "${outfile}" || true

                n=$(grep -c "^>" "${outfile}" || true)
                echo -e "${sp}\t${marker}\t${n}\t${outfile}" >> selected_other_refs/selection_log.tsv
            else
                : > "${outfile}"
                echo -e "${sp}\t${marker}\t0\tmissing" >> selected_other_refs/selection_log.tsv
            fi
        done
    done < other_species.list
}

merge_final_marker_refs() {
    log "Merging final marker references..."

    mkdir -p "${FINAL_MARKER_REF_DIR}"

    cat selected_marker_refs/TEF1_ref.selected.fasta selected_other_refs/TEF1/*.fasta \
        > "${FINAL_MARKER_REF_DIR}/TEF1_ref.final.fasta"

    cat selected_marker_refs/ACT_ref.selected.fasta selected_other_refs/ACT/*.fasta \
        > "${FINAL_MARKER_REF_DIR}/ACT_ref.final.fasta"

    cat selected_marker_refs/CAL_ref.selected.fasta selected_other_refs/CAL/*.fasta \
        > "${FINAL_MARKER_REF_DIR}/CAL_ref.final.fasta"

    cat selected_marker_refs/HIS3_ref.selected.fasta selected_other_refs/HIS3/*.fasta \
        > "${FINAL_MARKER_REF_DIR}/HIS3_ref.final.fasta"

    cat selected_marker_refs/ITS_ref.selected.fasta selected_other_refs/ITS/*.fasta \
        > "${FINAL_MARKER_REF_DIR}/ITS_ref.final.fasta"

    for marker in ${MARKERS}; do
        echo "===== ${marker} ====="
        grep -c "^>" "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" || true
        grep "^>" "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" || true
    done
}

main() {
    check_dependencies
    extract_headers
    build_candidate_tables
    filter_candidates
    prepare_core_accessions
    build_all_marker_pool
    extract_core_refs
    split_core_refs_by_marker
    select_other_species_refs
    merge_final_marker_refs

    log "Reference marker selection completed."
}

main "$@"
