# Fungal protein functional annotation pipeline

This repository contains scripts used for functional annotation of predicted proteins from two fungal genomes, referred to as CM and CZ.

## Overview

The pipeline performs the following analyses:

1. InterProScan domain annotation.
2. KOfamScan KEGG orthology annotation.
3. eggNOG-mapper functional annotation.
4. DIAMOND BLASTP against NCBI RefSeq fungi proteins.
5. dbCAN/CAZy-based carbohydrate-active enzyme and CWDE annotation.
6. MEROPS peptidase annotation.
7. PHI-base and UniProt fungi pathogenicity-related annotation.
8. Orthogroup-level integration of all functional annotation results.
9. Construction of a four-class CWDE table for downstream expression analysis.

## Input files

The following input files are required:

- `CM_proteins_cleaned.faa`
- `CZ_proteins_cleaned.faa`
- `CM_CZ_OG.expanded.csv` or `Complete_CM_CZ_Ortholog_Analysis.xlsx`

## Configuration

Edit `config.sh` before running the pipeline:

```bash
nano config.sh
```

Required database paths include:

- KOfam database
- eggNOG database
- InterProScan executable
- NCBI RefSeq fungi DIAMOND database
- dbCAN HMM database
- CAZy DIAMOND database
- MEROPS DIAMOND database
- PHI-base DIAMOND database
- UniProt fungi DIAMOND database

## Usage

Run individual steps:

```bash
bash 01_run_interpro_kofam_eggnog.sh -g both -s all -t 64
bash 02_run_ncbi_fungi_blast.sh
bash 03_run_dbcan_cwde.sh
bash 04_run_merops.sh
bash 05_run_phi_uniprot.sh
bash 07_run_integration.sh
python 08_build_cwde_4class_table.py
```

Or run all steps:

```bash
bash run_all_functional_annotation.sh
```

## Output directories

```text
annotation_results/
  interpro/
  kofam/
  eggnog/

orthofinder_input/
  blast_results/

dbcan_annotation/
merops_annotation/
phi_annotation/
integrated_annotation_results/
logs/
```

## Main outputs

Integrated annotation table:

```text
integrated_annotation_results/CM_CZ_OG_integrated_all_annotations.xlsx
integrated_annotation_results/CM_CZ_OG_integrated_all_annotations.csv
```

Specialized tables:

```text
integrated_annotation_results/CWDE_orthogroups.xlsx
integrated_annotation_results/Peptidase_orthogroups.xlsx
integrated_annotation_results/Pathogenicity_orthogroups.xlsx
```

CWDE four-class table:

```text
dbcan_annotation/CM_CZ_OG_expanded_CWDE_4class_primary_table.simple.xlsx
```

## Notes

Protein identifiers may appear in different formats across annotation tools, such as:

- `transcript:KAFxxxx`
- `transcript_KAFxxxx`
- `KAFxxxx`
- `KAMxxxx`

The integration script uses conservative alias matching to improve annotation recovery while avoiding inappropriate modification of CZ-specific KAM identifiers.
