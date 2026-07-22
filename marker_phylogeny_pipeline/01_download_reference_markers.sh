#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Download Cercospora/Pseudocercospora marker sequences from NCBI nucleotide
# using Entrez Direct.
#
# Markers:
#   ITS, TEF1, ACT, CAL, HIS3, TUB2
#
# Output:
#   selected_downloads/{ITS,TEF1,ACT,CAL,HIS3,TUB2}/*.fasta
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
    command -v esearch >/dev/null 2>&1 || die "esearch not found. Install Entrez Direct."
    command -v efetch >/dev/null 2>&1 || die "efetch not found. Install Entrez Direct."

    log "esearch: $(which esearch)"
    log "efetch : $(which efetch)"
}

prepare_species_list() {
    log "Preparing species list..."

    cat > "${SPECIES_LIST}" <<'LIST'
Cercospora zeae-maydis
Cercospora zeina
Cercospora beticola
Cercospora kikuchii
Cercospora sojina
Cercospora apii
Cercospora canescens
Cercospora nicotianae
Cercospora arachidicola
Cercospora coffeicola
Pseudocercospora fijiensis
LIST

    cat "${SPECIES_LIST}"
}

download_marker() {
    local marker="$1"
    local query_terms="$2"

    log "Downloading marker: ${marker}"

    mkdir -p "${SELECTED_DOWNLOADS_DIR}/${marker}"

    while read -r sp; do
        [[ -z "${sp}" ]] && continue

        local name
        name=$(echo "${sp}" | sed 's/ /_/g; s/\./_/g')

        local outfile="${SELECTED_DOWNLOADS_DIR}/${marker}/${name}_${marker}.fasta"

        log "Downloading ${marker} for ${sp}"

        esearch -db nucleotide \
            -query "\"${sp}\"[Organism] AND (${query_terms}) NOT genome NOT chromosome NOT scaffold" </dev/null \
        | efetch -format fasta > "${outfile}" || true

        local n
        n=$(grep -c "^>" "${outfile}" || true)

        log "Saved ${n} sequences to ${outfile}"

        sleep 1
    done < "${SPECIES_LIST}"
}

count_downloads() {
    log "Counting downloaded sequences..."

    for marker in TEF1 ACT CAL HIS3 TUB2 ITS; do
        echo "========== ${marker} =========="
        for f in "${SELECTED_DOWNLOADS_DIR}/${marker}"/*.fasta; do
            [[ -e "${f}" ]] || continue
            n=$(grep -c "^>" "${f}" || true)
            echo -e "${n}\t${f}"
        done
    done
}

main() {
    check_dependencies
    prepare_species_list

    mkdir -p "${SELECTED_DOWNLOADS_DIR}"/{TEF1,ACT,CAL,HIS3,TUB2,ITS}

    download_marker "TEF1" "\"translation elongation factor 1-alpha\" OR \"translation elongation factor 1 alpha\" OR \"elongation factor 1-alpha\" OR tef1 OR \"tef1-alpha\" OR \"EF1-alpha\""
    download_marker "ACT"  "actin OR ACT OR act1"
    download_marker "CAL"  "calmodulin OR CAL OR cmdA"
    download_marker "HIS3" "\"histone H3\" OR HIS3 OR his3"
    download_marker "TUB2" "\"beta-tubulin\" OR \"beta tubulin\" OR tub2 OR benA"
    download_marker "ITS"  "\"internal transcribed spacer\" OR ITS OR \"ribosomal RNA\""

    count_downloads

    log "Marker download completed."
}

main "$@"
