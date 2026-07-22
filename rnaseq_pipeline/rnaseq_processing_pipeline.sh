#!/usr/bin/env bash

# ==============================================================================
# RNA-seq preprocessing, mapping, host/pathogen BAM splitting, and read counting
#
# Steps:
#   1. Run FastQC on raw reads
#   2. Trim reads using Trimmomatic
#   3. Count reads after trimming
#   4. Build HISAT2 indexes
#   5. Map reads with HISAT2
#   6. Convert GFF3 to GTF
#   7. Prefix pathogen FASTA/GTF identifiers
#   8. Split combined-reference BAM files into host- and pathogen-derived BAMs
#   9. Quantify reads using featureCounts
#  10. Generate MultiQC reports
#
# Usage:
#   bash rnaseq_processing_pipeline.sh
#
# Requirements:
#   fastqc
#   multiqc
#   trimmomatic
#   hisat2
#   samtools
#   seqkit
#   gffread
#   featureCounts
#
# Author: Zhaoreng Li et al.
# ==============================================================================

set -euo pipefail

# Load configuration
CONFIG_FILE="./config.sh"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: config.sh not found."
    exit 1
fi

source "${CONFIG_FILE}"


# ==============================================================================
# Helper functions
# ==============================================================================

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

check_file() {
    local file="$1"
    if [[ ! -f "${file}" ]]; then
        echo "ERROR: Required file not found: ${file}"
        exit 1
    fi
}

check_dir() {
    local dir="$1"
    if [[ ! -d "${dir}" ]]; then
        echo "ERROR: Required directory not found: ${dir}"
        exit 1
    fi
}

make_dirs() {
    mkdir -p \
        "${FASTQC_PRETRIM_DIR}" \
        "${TRIMMED_DIR}" \
        "${FASTQC_POSTTRIM_DIR}" \
        "${GENOME_DIR}" \
        "${INDEX_DIR}" \
        "${BAM_DIR}" \
        "${MAPPING_STAT_DIR}" \
        "${SPLIT_BAM_DIR}/maize" \
        "${SPLIT_BAM_DIR}/cm" \
        "${SPLIT_BAM_DIR}/cz" \
        "${QUANT_DIR}" \
        "${LOG_DIR}"
}


# ==============================================================================
# Step 0. Check required input files
# ==============================================================================

check_inputs() {
    log_msg "Checking input files and directories..."

    check_dir "${RAW_DATA_DIR}"

    check_file "${MAIZE_FA}"
    check_file "${CM_FA}"
    check_file "${CZ_FA}"

    check_file "${MAIZE_GFF3}"
    check_file "${CM_GFF3}"
    check_file "${CZ_GFF3}"

    check_file "${ADAPTER_FA}"

    log_msg "Input check completed."
}


# ==============================================================================
# Step 1. FastQC before trimming
# ==============================================================================

run_fastqc_pretrim() {
    log_msg "Running FastQC on raw FASTQ files..."

    find "${RAW_DATA_DIR}" -type f -name "*.fastq.gz" | \
        xargs -I {} -P "${MAX_FASTQC_JOBS}" bash -c '
            fastqc -t "$0" --nogroup --kmers 7 -o "$1" "$2"
        ' "${THREADS_FASTQC}" "${FASTQC_PRETRIM_DIR}" {}

    multiqc "${FASTQC_PRETRIM_DIR}" \
        -n "${PRETRIM_MULTIQC_NAME}" \
        --force \
        -o "${FASTQC_PRETRIM_DIR}"

    log_msg "Pre-trim FastQC and MultiQC completed."
}


# ==============================================================================
# Step 2. Trim reads using Trimmomatic
# ==============================================================================

run_trimmomatic() {
    log_msg "Trimming reads using Trimmomatic..."

    for R1 in "${RAW_DATA_DIR}"/*_R1.fastq.gz; do
        [[ -e "${R1}" ]] || continue

        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

        if [[ ! -f "${R2}" ]]; then
            echo "WARNING: Missing paired file for ${R1}: ${R2}"
            continue
        fi

        BASE=$(basename "${R1}" _R1.fastq.gz)

        log_msg "Trimming sample: ${BASE}"

        trimmomatic PE \
            -threads "${THREADS_TRIMMOMATIC}" \
            -phred33 \
            "${R1}" "${R2}" \
            "${TRIMMED_DIR}/${BASE}_R1_paired.fastq.gz" \
            "${TRIMMED_DIR}/${BASE}_R1_unpaired.fastq.gz" \
            "${TRIMMED_DIR}/${BASE}_R2_paired.fastq.gz" \
            "${TRIMMED_DIR}/${BASE}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:"${ADAPTER_FA}":2:30:10 \
            SLIDINGWINDOW:4:15 \
            LEADING:10 \
            TRAILING:10 \
            MINLEN:50 \
            AVGQUAL:20 \
            > "${LOG_DIR}/${BASE}.trimmomatic.log" 2>&1
    done

    log_msg "Trimmomatic trimming completed."
}


# ==============================================================================
# Step 3. FastQC after trimming
# ==============================================================================

run_fastqc_posttrim() {
    log_msg "Running FastQC on trimmed paired FASTQ files..."

    find "${TRIMMED_DIR}" -type f -name "*_paired.fastq.gz" | \
        xargs -I {} -P "${MAX_FASTQC_JOBS}" bash -c '
            fastqc -t "$0" --nogroup --kmers 7 -o "$1" "$2"
        ' "${THREADS_FASTQC}" "${FASTQC_POSTTRIM_DIR}" {}

    multiqc "${FASTQC_POSTTRIM_DIR}" \
        -n "${POSTTRIM_MULTIQC_NAME}" \
        --force \
        -o "${FASTQC_POSTTRIM_DIR}"

    log_msg "Post-trim FastQC and MultiQC completed."
}


# ==============================================================================
# Step 4. Count reads after trimming
# ==============================================================================

count_trimmed_reads() {
    log_msg "Counting reads after trimming..."

    local outfile="${TRIMMED_DIR}/read_count_summary.txt"

    echo -e "Sample\tR1_reads\tR2_reads\tTotal_reads" > "${outfile}"

    for R1 in "${TRIMMED_DIR}"/*_R1_paired.fastq.gz; do
        [[ -e "${R1}" ]] || continue

        R2="${R1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}"
        SAMPLE=$(basename "${R1}" _R1_paired.fastq.gz)

        if [[ ! -f "${R2}" ]]; then
            echo "WARNING: Missing paired file for ${SAMPLE}: ${R2}"
            continue
        fi

        R1_LINES=$(zcat "${R1}" | wc -l)
        R2_LINES=$(zcat "${R2}" | wc -l)

        R1_READS=$((R1_LINES / 4))
        R2_READS=$((R2_LINES / 4))
        TOTAL=$((R1_READS + R2_READS))

        echo -e "${SAMPLE}\t${R1_READS}\t${R2_READS}\t${TOTAL}" >> "${outfile}"
    done

    log_msg "Trimmed read count table saved: ${outfile}"
}


# ==============================================================================
# Step 5. Prepare prefixed pathogen references and combined references
# ==============================================================================

prepare_references() {
    log_msg "Preparing prefixed pathogen references and combined host-pathogen references..."

    # Add prefixes to pathogen FASTA sequence IDs
    seqkit replace -p "^(.+)$" -r "CM_\$1" "${CM_FA}" > "${CM_PREFIXED_FA}"
    seqkit replace -p "^(.+)$" -r "CZ_\$1" "${CZ_FA}" > "${CZ_PREFIXED_FA}"

    # Concatenate maize and pathogen references
    cat "${MAIZE_FA}" "${CM_PREFIXED_FA}" > "${ZM_CM_FA}"
    cat "${MAIZE_FA}" "${CZ_PREFIXED_FA}" > "${ZM_CZ_FA}"

    log_msg "Reference preparation completed."
}


# ==============================================================================
# Step 6. Build HISAT2 indexes
# ==============================================================================

build_hisat2_indexes() {
    log_msg "Building HISAT2 indexes..."

    hisat2-build -p "${THREADS_HISAT2_BUILD}" "${MAIZE_FA}" "${INDEX_ZM}"
    hisat2-build -p "${THREADS_HISAT2_BUILD}" "${CM_PREFIXED_FA}" "${INDEX_CM}"
    hisat2-build -p "${THREADS_HISAT2_BUILD}" "${CZ_PREFIXED_FA}" "${INDEX_CZ}"
    hisat2-build -p "${THREADS_HISAT2_BUILD}" "${ZM_CM_FA}" "${INDEX_ZM_CM}"
    hisat2-build -p "${THREADS_HISAT2_BUILD}" "${ZM_CZ_FA}" "${INDEX_ZM_CZ}"

    log_msg "HISAT2 index construction completed."
}


# ==============================================================================
# Step 7. Determine HISAT2 index according to sample name
# ==============================================================================

select_index_for_sample() {
    local sample="$1"

    # Maize-only control samples
    if [[ "${sample}" == CK-* || "${sample}" == CK_* ]]; then
        echo "${INDEX_ZM}"

    # In vitro CM samples
    elif [[ "${sample}" == CM_* ]]; then
        echo "${INDEX_CM}"

    # In vitro CZ samples
    elif [[ "${sample}" == CZ_* ]]; then
        echo "${INDEX_CZ}"

    # CM-inoculated maize samples
    elif [[ "${sample}" == CM0_* || "${sample}" == CM1-* || "${sample}" == CM2-* || "${sample}" == CM3-* || \
            "${sample}" == CM1_* || "${sample}" == CM2_* || "${sample}" == CM3_* ]]; then
        echo "${INDEX_ZM_CM}"

    # CZ-inoculated maize samples
    elif [[ "${sample}" == CZ0_* || "${sample}" == CZ1-* || "${sample}" == CZ2-* || "${sample}" == CZ3-* || \
            "${sample}" == CZ1_* || "${sample}" == CZ2_* || "${sample}" == CZ3_* ]]; then
        echo "${INDEX_ZM_CZ}"

    else
        echo "UNKNOWN"
    fi
}


# ==============================================================================
# Step 8. Map reads using HISAT2
# ==============================================================================

run_hisat2_mapping() {
    log_msg "Mapping trimmed reads with HISAT2..."

    for R1 in "${TRIMMED_DIR}"/*_R1_paired.fastq.gz; do
        [[ -e "${R1}" ]] || continue

        R2="${R1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}"
        SAMPLE=$(basename "${R1}" _R1_paired.fastq.gz)

        if [[ ! -f "${R2}" ]]; then
            echo "WARNING: Missing paired file for ${SAMPLE}: ${R2}"
            continue
        fi

        INDEX=$(select_index_for_sample "${SAMPLE}")

        if [[ "${INDEX}" == "UNKNOWN" ]]; then
            echo "WARNING: Unknown sample type: ${SAMPLE}. Skipping."
            continue
        fi

        log_msg "Mapping sample ${SAMPLE} using index ${INDEX}"

        hisat2 \
            -p "${THREADS_HISAT2}" \
            -x "${INDEX}" \
            -1 "${R1}" \
            -2 "${R2}" \
            2> "${MAPPING_STAT_DIR}/${SAMPLE}.hisat2.log" | \
            samtools view -@ "${THREADS_SAMTOOLS}" -bS - | \
            samtools sort -@ "${THREADS_SAMTOOLS}" \
                -o "${BAM_DIR}/${SAMPLE}.sorted.bam"

        samtools index "${BAM_DIR}/${SAMPLE}.sorted.bam"
    done

    log_msg "HISAT2 mapping completed."
}


# ==============================================================================
# Step 9. Convert GFF3 annotation to GTF
# ==============================================================================

convert_gff3_to_gtf() {
    log_msg "Converting GFF3 files to GTF format..."

    gffread "${MAIZE_GFF3}" -T -o "${MAIZE_GTF}"
    gffread "${CM_GFF3}" -T -o "${CM_GTF}"
    gffread "${CZ_GFF3}" -T -o "${CZ_GTF}"

    log_msg "GFF3 to GTF conversion completed."
}


# ==============================================================================
# Step 10. Prefix pathogen GTF reference sequence names
# ==============================================================================

prefix_pathogen_gtf() {
    log_msg "Adding prefixes to pathogen GTF reference sequence names..."

    awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$1="CM_"$1; print}' \
        "${CM_GTF}" > "${CM_PREFIXED_GTF}"

    awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$1="CZ_"$1; print}' \
        "${CZ_GTF}" > "${CZ_PREFIXED_GTF}"

    log_msg "Prefixed pathogen GTF files generated:"
    log_msg "  ${CM_PREFIXED_GTF}"
    log_msg "  ${CZ_PREFIXED_GTF}"
}


# ==============================================================================
# Step 11. Generate reference sequence lists for BAM splitting
# ==============================================================================

generate_ref_lists_for_splitting() {
    log_msg "Generating reference sequence lists for BAM splitting..."

    local cm_example_bam
    local cz_example_bam

    cm_example_bam=$(find "${BAM_DIR}" -name "CM0_*.sorted.bam" -o -name "CM0-*.sorted.bam" | head -n 1 || true)
    cz_example_bam=$(find "${BAM_DIR}" -name "CZ0_*.sorted.bam" -o -name "CZ0-*.sorted.bam" | head -n 1 || true)

    if [[ -z "${cm_example_bam}" ]]; then
        echo "ERROR: No CM combined-reference BAM found for generating reference lists."
        exit 1
    fi

    if [[ -z "${cz_example_bam}" ]]; then
        echo "ERROR: No CZ combined-reference BAM found for generating reference lists."
        exit 1
    fi

    samtools idxstats "${cm_example_bam}" | cut -f1 | grep -v '^\*$' | grep '^CM_' > "${CM_REF_LIST}"
    samtools idxstats "${cm_example_bam}" | cut -f1 | grep -v '^\*$' | grep -v '^CM_' > "${MAIZE_FROM_CM_REF_LIST}"

    samtools idxstats "${cz_example_bam}" | cut -f1 | grep -v '^\*$' | grep '^CZ_' > "${CZ_REF_LIST}"
    samtools idxstats "${cz_example_bam}" | cut -f1 | grep -v '^\*$' | grep -v '^CZ_' > "${MAIZE_FROM_CZ_REF_LIST}"

    log_msg "Reference lists generated:"
    log_msg "  ${CM_REF_LIST}"
    log_msg "  ${CZ_REF_LIST}"
    log_msg "  ${MAIZE_FROM_CM_REF_LIST}"
    log_msg "  ${MAIZE_FROM_CZ_REF_LIST}"
}


# ==============================================================================
# Step 12. Split BAM files by host/pathogen reference sequences
# ==============================================================================

extract_by_list() {
    local bam="$1"
    local list="$2"
    local outbam="$3"

    if [[ ! -s "${list}" ]]; then
        echo "ERROR: Reference list does not exist or is empty: ${list}"
        exit 1
    fi

    mapfile -t refs < "${list}"

    samtools view -b "${bam}" "${refs[@]}" > "${outbam}"
    samtools index "${outbam}"
}


split_bam_files() {
    log_msg "Splitting BAM files into maize- and pathogen-derived BAMs..."

    for BAM in "${BAM_DIR}"/*.sorted.bam; do
        [[ -e "${BAM}" ]] || continue

        SAMPLE=$(basename "${BAM}" .sorted.bam)

        log_msg "Processing BAM: ${SAMPLE}"

        if [[ "${SAMPLE}" == CK-* || "${SAMPLE}" == CK_* ]]; then
            cp "${BAM}" "${SPLIT_BAM_DIR}/maize/${SAMPLE}.maize.bam"
            samtools index "${SPLIT_BAM_DIR}/maize/${SAMPLE}.maize.bam"

        elif [[ "${SAMPLE}" == CM_* ]]; then
            cp "${BAM}" "${SPLIT_BAM_DIR}/cm/${SAMPLE}.cm.bam"
            samtools index "${SPLIT_BAM_DIR}/cm/${SAMPLE}.cm.bam"

        elif [[ "${SAMPLE}" == CZ_* ]]; then
            cp "${BAM}" "${SPLIT_BAM_DIR}/cz/${SAMPLE}.cz.bam"
            samtools index "${SPLIT_BAM_DIR}/cz/${SAMPLE}.cz.bam"

        elif [[ "${SAMPLE}" == CM0_* || "${SAMPLE}" == CM1-* || "${SAMPLE}" == CM2-* || "${SAMPLE}" == CM3-* || \
                "${SAMPLE}" == CM1_* || "${SAMPLE}" == CM2_* || "${SAMPLE}" == CM3_* ]]; then
            extract_by_list "${BAM}" "${MAIZE_FROM_CM_REF_LIST}" "${SPLIT_BAM_DIR}/maize/${SAMPLE}.maize.bam"
            extract_by_list "${BAM}" "${CM_REF_LIST}" "${SPLIT_BAM_DIR}/cm/${SAMPLE}.cm.bam"

        elif [[ "${SAMPLE}" == CZ0_* || "${SAMPLE}" == CZ1-* || "${SAMPLE}" == CZ2-* || "${SAMPLE}" == CZ3-* || \
                "${SAMPLE}" == CZ1_* || "${SAMPLE}" == CZ2_* || "${SAMPLE}" == CZ3_* ]]; then
            extract_by_list "${BAM}" "${MAIZE_FROM_CZ_REF_LIST}" "${SPLIT_BAM_DIR}/maize/${SAMPLE}.maize.bam"
            extract_by_list "${BAM}" "${CZ_REF_LIST}" "${SPLIT_BAM_DIR}/cz/${SAMPLE}.cz.bam"

        else
            echo "WARNING: Unknown sample type during BAM splitting: ${SAMPLE}"
        fi
    done

    log_msg "BAM splitting completed."
}


# ==============================================================================
# Step 13. Quantify reads using featureCounts
# ==============================================================================

run_featurecounts() {
    log_msg "Running featureCounts..."

    # Maize gene-level counting
    featureCounts \
        -T "${THREADS_FEATURECOUNTS}" \
        -p \
        -t exon \
        -g gene_id \
        -a "${MAIZE_GTF}" \
        -o "${MAIZE_COUNTS}" \
        "${SPLIT_BAM_DIR}"/maize/*.maize.bam

    # CM pathogen counting
    featureCounts \
        -T "${THREADS_FEATURECOUNTS}" \
        -p \
        -t exon \
        -g transcript_id \
        -a "${CM_PREFIXED_GTF}" \
        -o "${CM_COUNTS}" \
        "${SPLIT_BAM_DIR}"/cm/*.cm.bam

    # CZ pathogen counting
    featureCounts \
        -T "${THREADS_FEATURECOUNTS}" \
        -p \
        -t exon \
        -g transcript_id \
        -a "${CZ_PREFIXED_GTF}" \
        -o "${CZ_COUNTS}" \
        "${SPLIT_BAM_DIR}"/cz/*.cz.bam

    log_msg "featureCounts quantification completed."
}


# ==============================================================================
# Step 14. Run MultiQC on logs
# ==============================================================================

run_multiqc_all() {
    log_msg "Running MultiQC on the whole project..."

    multiqc "${PROJECT_DIR}" \
        -n "rnaseq_processing_multiqc_report" \
        --force \
        -o "${PROJECT_DIR}"

    log_msg "Final MultiQC report completed."
}


# ==============================================================================
# Main workflow
# ==============================================================================

main() {
    log_msg "RNA-seq processing pipeline started."

    make_dirs
    check_inputs

    run_fastqc_pretrim
    run_trimmomatic
    run_fastqc_posttrim
    count_trimmed_reads

    prepare_references
    build_hisat2_indexes

    run_hisat2_mapping

    convert_gff3_to_gtf
    prefix_pathogen_gtf
    generate_ref_lists_for_splitting
    split_bam_files

    run_featurecounts
    run_multiqc_all

    log_msg "RNA-seq processing pipeline completed successfully."
}

main "$@"
