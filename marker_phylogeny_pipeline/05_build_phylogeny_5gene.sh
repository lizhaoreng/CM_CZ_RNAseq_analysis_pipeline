#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Build multilocus phylogeny using ITS, TEF1, ACT, CAL, and HIS3.
#
# Input:
#   final_marker_refs/{ITS,TEF1,ACT,CAL,HIS3}_ref.final.fasta
#   own_marker_regions/renamed/{CM,CZ}_{TEF1,ACT,CAL,HIS3}.fasta
#   genome_marker_regions/renamed/{CM,CZ}_genome_{ITS,TEF1,ACT,CAL,HIS3}.fasta
#
# Output:
#   final_marker_with_own/
#   alignments/
#   alignments_trimmed/
#   concatenated/
#   iqtree_out/
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
    for cmd in seqkit mafft trimal AMAS.py iqtree; do
        command -v "$cmd" >/dev/null 2>&1 || die "$cmd not found"
    done
}

check_inputs() {
    for marker in ${MARKERS_ALL}; do
        ref="${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta"
        [[ -s "${ref}" ]] || die "Missing marker reference: ${ref}"
    done

    for marker in ${MARKERS_CODING}; do
        for sample in CM CZ; do
            f="${OWN_MARKER_DIR}/renamed/${sample}_${marker}.fasta"
            [[ -s "${f}" ]] || die "Missing Trinity marker file: ${f}"
        done
    done

    for marker in ${MARKERS_ALL}; do
        for sample in CM CZ; do
            f="${GENOME_MARKER_DIR}/renamed/${sample}_genome_${marker}.fasta"
            [[ -s "${f}" ]] || die "Missing genome marker file: ${f}"
        done
    done
}

merge_markers() {
    log "Merging reference, Trinity-derived, and genome-derived markers..."

    mkdir -p "${FINAL_MARKER_WITH_OWN_DIR}"

    cat "${FINAL_MARKER_REF_DIR}/ITS_ref.final.fasta" \
        "${GENOME_MARKER_DIR}/renamed/CM_genome_ITS.fasta" \
        "${GENOME_MARKER_DIR}/renamed/CZ_genome_ITS.fasta" \
        > "${FINAL_MARKER_WITH_OWN_DIR}/ITS_all.fasta"

    for marker in ${MARKERS_CODING}; do
        cat "${FINAL_MARKER_REF_DIR}/${marker}_ref.final.fasta" \
            "${OWN_MARKER_DIR}/renamed/CM_${marker}.fasta" \
            "${OWN_MARKER_DIR}/renamed/CZ_${marker}.fasta" \
            "${GENOME_MARKER_DIR}/renamed/CM_genome_${marker}.fasta" \
            "${GENOME_MARKER_DIR}/renamed/CZ_genome_${marker}.fasta" \
            > "${FINAL_MARKER_WITH_OWN_DIR}/${marker}_all.fasta"
    done

    for marker in ${MARKERS_ALL}; do
        echo "===== ${marker} all ====="
        grep -c "^>" "${FINAL_MARKER_WITH_OWN_DIR}/${marker}_all.fasta" || true
    done
}

rename_by_taxon() {
    log "Renaming sequences by taxon..."

    mkdir -p "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon"

    # ==========================================================================
    # ITS
    # ==========================================================================
    seqkit replace -p "^OR945717.1.*" -r "CM_Jilin_this_study" "${FINAL_MARKER_WITH_OWN_DIR}/ITS_all.fasta" \
    | seqkit replace -p "^PV133732.1.*" -r "CZ_Jilin_this_study" \
    | seqkit replace -p "^CM_genome_ITS.*" -r "C_zeae_maydis_SCOH1_5_genome" \
    | seqkit replace -p "^CZ_genome_ITS.*" -r "C_zeina_CMW25467_genome" \
    | seqkit replace -p "^KM087697.1.*" -r "C_zeae_maydis_2197" \
    | seqkit replace -p "^DQ185080.1.*" -r "C_zeae_maydis_CBS117763" \
    | seqkit replace -p "^NR_111205.1.*" -r "C_zeina_CPC11995_extype" \
    | seqkit replace -p "^EU569227.1.*" -r "C_zeina_CMW25467" \
    | seqkit replace -p "^PZ556481.1.*" -r "C_apii" \
    | seqkit replace -p "^PX898671.1.*" -r "C_asparagi" \
    | seqkit replace -p "^MH855153.1.*" -r "C_beticola" \
    | seqkit replace -p "^PZ322670.1.*" -r "C_cf_flagellaris" \
    | seqkit replace -p "^PX754784.1.*" -r "C_flagellaris" \
    | seqkit replace -p "^MH854904.1.*" -r "C_kikuchii" \
    | seqkit replace -p "^PX923175.1.*" -r "C_nicotianae" \
    | seqkit replace -p "^PV089366.1.*" -r "C_sojina" \
    | seqkit replace -p "^PZ017509.1.*" -r "P_fijiensis" \
    | seqkit replace -p "^PX945723.1.*" -r "P_griseola" \
    > "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/ITS.fasta"

    # ==========================================================================
    # TEF1
    # ==========================================================================
    seqkit replace -p "^KM087649.1.*" -r "C_zeae_maydis_2197" "${FINAL_MARKER_WITH_OWN_DIR}/TEF1_all.fasta" \
    | seqkit replace -p "^DQ185092.1.*" -r "C_zeae_maydis_CBS117763" \
    | seqkit replace -p "^EU569218.1.*" -r "C_zeina_CMW25467" \
    | seqkit replace -p "^DQ185094.1.*" -r "C_zeina_CPC11998" \
    | seqkit replace -p "^DQ185093.1.*" -r "C_zeina_CPC11995_extype" \
    | seqkit replace -p "^KF253244.1.*" -r "C_apii" \
    | seqkit replace -p "^PX906524.1.*" -r "C_asparagi" \
    | seqkit replace -p "^PX104760.1.*" -r "C_beticola" \
    | seqkit replace -p "^PV335258.1.*" -r "C_cf_flagellaris" \
    | seqkit replace -p "^PX072464.1.*" -r "C_flagellaris" \
    | seqkit replace -p "^PX515994.1.*" -r "C_kikuchii" \
    | seqkit replace -p "^MT013706.1.*" -r "C_nicotianae" \
    | seqkit replace -p "^OP594329.1.*" -r "C_sojina" \
    | seqkit replace -p "^PP404794.1.*" -r "P_fijiensis" \
    | seqkit replace -p "^PX985209.1.*" -r "P_griseola" \
    | seqkit replace -p "^CM_TEF1.*" -r "CM_Jilin_this_study" \
    | seqkit replace -p "^CZ_TEF1.*" -r "CZ_Jilin_this_study" \
    | seqkit replace -p "^CM_genome_TEF1.*" -r "C_zeae_maydis_SCOH1_5_genome" \
    | seqkit replace -p "^CZ_genome_TEF1.*" -r "C_zeina_CMW25467_genome" \
    > "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/TEF1.fasta"

    # ==========================================================================
    # ACT
    # ==========================================================================
    seqkit replace -p "^DQ185104.1.*" -r "C_zeae_maydis_CBS117763" "${FINAL_MARKER_WITH_OWN_DIR}/ACT_all.fasta" \
    | seqkit replace -p "^DQ185106.1.*" -r "C_zeina_CPC11998" \
    | seqkit replace -p "^DQ185105.1.*" -r "C_zeina_CPC11995_extype" \
    | seqkit replace -p "^PX663050.1.*" -r "C_apii" \
    | seqkit replace -p "^PX906515.1.*" -r "C_asparagi" \
    | seqkit replace -p "^PX794519.1.*" -r "C_beticola" \
    | seqkit replace -p "^OR548247.1.*" -r "C_cf_flagellaris" \
    | seqkit replace -p "^PX072446.1.*" -r "C_flagellaris" \
    | seqkit replace -p "^PP918308.1.*" -r "C_kikuchii" \
    | seqkit replace -p "^MT013583.1.*" -r "C_nicotianae" \
    | seqkit replace -p "^MZ456945.1.*" -r "C_sojina" \
    | seqkit replace -p "^PX210785.1.*" -r "P_fijiensis" \
    | seqkit replace -p "^PX985207.1.*" -r "P_griseola" \
    | seqkit replace -p "^CM_ACT.*" -r "CM_Jilin_this_study" \
    | seqkit replace -p "^CZ_ACT.*" -r "CZ_Jilin_this_study" \
    | seqkit replace -p "^CM_genome_ACT.*" -r "C_zeae_maydis_SCOH1_5_genome" \
    | seqkit replace -p "^CZ_genome_ACT.*" -r "C_zeina_CMW25467_genome" \
    > "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/ACT.fasta"

    # ==========================================================================
    # CAL
    # ==========================================================================
    seqkit replace -p "^OQ773855.1.*" -r "C_zeae_maydis_Czm_20_107" "${FINAL_MARKER_WITH_OWN_DIR}/CAL_all.fasta" \
    | seqkit replace -p "^DQ185118.1.*" -r "C_zeina_CPC11998" \
    | seqkit replace -p "^DQ185117.1.*" -r "C_zeina_CPC11995_extype" \
    | seqkit replace -p "^PX663051.1.*" -r "C_apii" \
    | seqkit replace -p "^PX906518.1.*" -r "C_asparagi" \
    | seqkit replace -p "^PX104878.1.*" -r "C_beticola" \
    | seqkit replace -p "^PV335268.1.*" -r "C_cf_flagellaris" \
    | seqkit replace -p "^PX072452.1.*" -r "C_flagellaris" \
    | seqkit replace -p "^PX515993.1.*" -r "C_kikuchii" \
    | seqkit replace -p "^OQ190497.1.*" -r "C_nicotianae" \
    | seqkit replace -p "^OP594347.1.*" -r "C_sojina" \
    | seqkit replace -p "^CM_CAL.*" -r "CM_Jilin_this_study" \
    | seqkit replace -p "^CZ_CAL.*" -r "CZ_Jilin_this_study" \
    | seqkit replace -p "^CM_genome_CAL.*" -r "C_zeae_maydis_SCOH1_5_genome" \
    | seqkit replace -p "^CZ_genome_CAL.*" -r "C_zeina_CMW25467_genome" \
    > "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/CAL.fasta"

    # ==========================================================================
    # HIS3
    # ==========================================================================
    seqkit replace -p "^KM087673.1.*" -r "C_zeae_maydis_2197" "${FINAL_MARKER_WITH_OWN_DIR}/HIS3_all.fasta" \
    | seqkit replace -p "^DQ185128.1.*" -r "C_zeae_maydis_CBS117763" \
    | seqkit replace -p "^DQ185130.1.*" -r "C_zeina_CPC11998" \
    | seqkit replace -p "^DQ185129.1.*" -r "C_zeina_CPC11995_extype" \
    | seqkit replace -p "^OQ790157.1.*" -r "C_apii" \
    | seqkit replace -p "^PX906521.1.*" -r "C_asparagi" \
    | seqkit replace -p "^PX104819.1.*" -r "C_beticola" \
    | seqkit replace -p "^PX091955.1.*" -r "C_cf_flagellaris" \
    | seqkit replace -p "^PX072458.1.*" -r "C_flagellaris" \
    | seqkit replace -p "^PX515995.1.*" -r "C_kikuchii" \
    | seqkit replace -p "^OQ241186.1.*" -r "C_nicotianae" \
    | seqkit replace -p "^MZ456955.1.*" -r "C_sojina" \
    | seqkit replace -p "^EU514370.1.*" -r "P_fijiensis" \
    | seqkit replace -p "^CM_HIS3.*" -r "CM_Jilin_this_study" \
    | seqkit replace -p "^CZ_HIS3.*" -r "CZ_Jilin_this_study" \
    | seqkit replace -p "^CM_genome_HIS3.*" -r "C_zeae_maydis_SCOH1_5_genome" \
    | seqkit replace -p "^CZ_genome_HIS3.*" -r "C_zeina_CMW25467_genome" \
    > "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/HIS3.fasta"
}

check_renamed_files() {
    log "Checking renamed marker files..."

    for marker in ${MARKERS_ALL}; do
        echo "===== ${marker} ====="
        grep -c "^>" "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/${marker}.fasta" || true
        grep "^>" "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/${marker}.fasta" || true
    done

    log "Checking possible unrenamed IDs..."

    for marker in ${MARKERS_ALL}; do
        echo "===== ${marker} possible unrenamed IDs ====="
        grep "^>" "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/${marker}.fasta" \
            | grep -E "^>[A-Z]{1,3}[0-9]|^>NR_|^>OQ|^>PX|^>DQ|^>KM|^>EU|^>MT|^>OP|^>PP|^>PV|^>MH|^>MZ|^>OR|^>KF" || true
    done

    log "Checking duplicate IDs..."

    for marker in ${MARKERS_ALL}; do
        echo "===== ${marker} duplicate IDs ====="
        grep "^>" "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/${marker}.fasta" \
            | sed 's/^>//' \
            | sort \
            | uniq -d || true
    done
}

run_alignment_and_trimming() {
    log "Running MAFFT alignment and trimAl trimming..."

    mkdir -p "${ALIGNMENT_DIR}" "${TRIMMED_ALIGNMENT_DIR}"

    for marker in ${MARKERS_ALL}; do
        mafft --auto --thread "${THREADS}" \
            "${FINAL_MARKER_WITH_OWN_DIR}/renamed_by_taxon/${marker}.fasta" \
            > "${ALIGNMENT_DIR}/${marker}.aln.fasta"

        trimal \
            -in "${ALIGNMENT_DIR}/${marker}.aln.fasta" \
            -out "${TRIMMED_ALIGNMENT_DIR}/${marker}.trimmed.fasta" \
            -automated1

        seqkit stats "${TRIMMED_ALIGNMENT_DIR}/${marker}.trimmed.fasta"
    done
}

concatenate_and_run_iqtree_5gene() {
    log "Concatenating five-gene matrix and running IQ-TREE..."

    mkdir -p "${CONCAT_DIR}" "${IQTREE_DIR}"

    AMAS.py concat \
        -i "${TRIMMED_ALIGNMENT_DIR}/ITS.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/TEF1.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/ACT.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/CAL.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/HIS3.trimmed.fasta" \
        -f fasta \
        -d dna \
        -p "${CONCAT_DIR}/partitions.txt" \
        -t "${CONCAT_DIR}/Cercospora_ITS_TEF1_ACT_CAL_HIS3_concat.fasta"

    cp "${CONCAT_DIR}/partitions.txt" "${CONCAT_DIR}/partitions_iqtree.txt"

    seqkit stats "${CONCAT_DIR}/Cercospora_ITS_TEF1_ACT_CAL_HIS3_concat.fasta"
    cat "${CONCAT_DIR}/partitions_iqtree.txt"

    iqtree \
        -s "${CONCAT_DIR}/Cercospora_ITS_TEF1_ACT_CAL_HIS3_concat.fasta" \
        -p "${CONCAT_DIR}/partitions_iqtree.txt" \
        -m MFP+MERGE \
        -B 1000 \
        -alrt 1000 \
        -T "${IQTREE_THREADS}" \
        -redo \
        --prefix "${IQTREE_DIR}/Cercospora_5gene_with_genomes"
}

concatenate_and_run_iqtree_noCAL() {
    log "Optional: building four-locus tree excluding CAL..."

    mkdir -p "${CONCAT_NOCAL_DIR}" "${IQTREE_NOCAL_DIR}"

    AMAS.py concat \
        -i "${TRIMMED_ALIGNMENT_DIR}/ITS.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/TEF1.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/ACT.trimmed.fasta" \
           "${TRIMMED_ALIGNMENT_DIR}/HIS3.trimmed.fasta" \
        -f fasta \
        -d dna \
        -p "${CONCAT_NOCAL_DIR}/partitions.txt" \
        -t "${CONCAT_NOCAL_DIR}/Cercospora_ITS_TEF1_ACT_HIS3_concat.fasta"

    cp "${CONCAT_NOCAL_DIR}/partitions.txt" "${CONCAT_NOCAL_DIR}/partitions_iqtree.txt"

    seqkit stats "${CONCAT_NOCAL_DIR}/Cercospora_ITS_TEF1_ACT_HIS3_concat.fasta"
    cat "${CONCAT_NOCAL_DIR}/partitions_iqtree.txt"

    iqtree \
        -s "${CONCAT_NOCAL_DIR}/Cercospora_ITS_TEF1_ACT_HIS3_concat.fasta" \
        -p "${CONCAT_NOCAL_DIR}/partitions_iqtree.txt" \
        -m MFP+MERGE \
        -B 1000 \
        -alrt 1000 \
        -T "${IQTREE_THREADS}" \
        -redo \
        --prefix "${IQTREE_NOCAL_DIR}/Cercospora_4locus_noCAL_with_genomes"
}

main() {
    check_dependencies
    check_inputs
    merge_markers
    rename_by_taxon
    check_renamed_files
    run_alignment_and_trimming
    concatenate_and_run_iqtree_5gene
    concatenate_and_run_iqtree_noCAL

    log "Phylogenetic analysis completed."
    log "Five-gene tree: ${IQTREE_DIR}/Cercospora_5gene_with_genomes.treefile"
    log "No-CAL tree: ${IQTREE_NOCAL_DIR}/Cercospora_4locus_noCAL_with_genomes.treefile"
}

main "$@"
