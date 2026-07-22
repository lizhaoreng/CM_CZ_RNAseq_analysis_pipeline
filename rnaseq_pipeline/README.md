# RNA-seq preprocessing and quantification pipeline

This repository contains the shell scripts used for RNA-seq quality control, trimming, read alignment, host/pathogen BAM splitting, and read counting in the comparative transcriptomic analysis of representative *Cercospora zeae-maydis* and *Cercospora zeina* isolates during maize infection.

## Overview

The pipeline performs the following steps:

1. Raw read quality control with FastQC and MultiQC.
2. Adapter and quality trimming using Trimmomatic.
3. Post-trimming quality control.
4. Read count summary for trimmed paired-end reads.
5. Construction of HISAT2 indexes for maize, pathogen, and combined host-pathogen references.
6. HISAT2 alignment of RNA-seq reads.
7. Conversion of GFF3 annotations to GTF format.
8. Addition of pathogen-specific prefixes to reference sequence IDs and GTF annotations.
9. Splitting combined-reference BAM files into maize-derived and pathogen-derived BAM files.
10. Gene or transcript-level read counting using featureCounts.
11. Final MultiQC report generation.

## Required software

The following programs are required:

- FastQC
- MultiQC
- Trimmomatic
- HISAT2
- SAMtools
- SeqKit
- gffread
- featureCounts/subread

A conda environment can be created using:

```bash
conda env create -f environment.yml
conda activate rnaseq_processing
```

## Input files

The following files should be specified in `config.sh`:

- Maize reference genome FASTA
- *C. zeae-maydis* reference genome FASTA
- *C. zeina* reference genome FASTA
- Maize GFF3 annotation
- *C. zeae-maydis* GFF3 annotation
- *C. zeina* GFF3 annotation
- Paired-end FASTQ files
- Trimmomatic adapter file

## FASTQ naming convention

The script assumes paired-end FASTQ files are named as:

```text
Sample_R1.fastq.gz
Sample_R2.fastq.gz
```

Example:

```text
CM0_1_R1.fastq.gz
CM0_1_R2.fastq.gz
CZ3_2_R1.fastq.gz
CZ3_2_R2.fastq.gz
```

## Sample type recognition

The script assigns references based on sample names:

| Sample pattern | Reference used |
|---|---|
| `CK-*` or `CK_*` | maize reference |
| `CM_*` | *C. zeae-maydis* pathogen reference |
| `CZ_*` | *C. zeina* pathogen reference |
| `CM0_*`, `CM1-*`, `CM2-*`, `CM3-*` | maize + *C. zeae-maydis* combined reference |
| `CZ0_*`, `CZ1-*`, `CZ2-*`, `CZ3-*` | maize + *C. zeina* combined reference |

Users should modify the `select_index_for_sample()` function if their sample naming scheme differs.

## Usage

Edit `config.sh` first:

```bash
nano config.sh
```

Then run:

```bash
bash rnaseq_processing_pipeline.sh
```

## Outputs

Main output directories include:

```text
fastqc_pretrim/
trimmed/
fastqc_posttrim/
genome/index/
mapping/bam/
mapping/stat/
mapping/split_bam/
featureCounts/
logs/
```

Main count tables:

```text
featureCounts/maize_counts.txt
featureCounts/cm_counts.txt
featureCounts/cz_counts.txt
```

## Notes

For host-pathogen mixed RNA-seq libraries, reads were aligned to combined maize-pathogen references. BAM files were then separated into maize-derived and pathogen-derived BAM files based on reference sequence prefixes. Pathogen genome sequence IDs were prefixed with `CM_` or `CZ_` before combined-reference construction to avoid sequence name conflicts.
