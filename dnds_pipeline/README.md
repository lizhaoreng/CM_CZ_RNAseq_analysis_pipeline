# dN/dS evolutionary analysis pipeline

This repository contains scripts used to estimate pairwise dN/dS ratios for single-copy orthologous genes and to compare evolutionary rates among candidate effectors, CWDEs, and background genes.

## Pipeline overview

The pipeline performs:

1. Ortholog table processing and FASTA identifier normalization.
2. Renaming of CDS and protein sequences.
3. Extraction of one-to-one orthologous CDS and protein sequence pairs.
4. Protein alignment using MAFFT.
5. Codon alignment using PAL2NAL.
6. Pairwise dN/dS estimation using PAML codeml.
7. Alignment quality assessment.
8. codeml result parsing and quality filtering.
9. Statistical comparison of evolutionary rates among functional gene categories.

## Requirements

The pipeline requires:

- Python >= 3.9
- pandas
- Biopython
- PyYAML
- MAFFT
- PAL2NAL
- PAML/codeml
- R >= 4.0
- tidyverse
- data.table
- ggpubr
- patchwork
- rstatix
- broom
- Cairo

A conda environment can be created using:

```bash
conda env create -f environment.yml
conda activate dnds_pipeline
```

## Input files

Expected input files are specified in `config.yaml`.

Main inputs:

```text
input/One_to_One_Orthologs.csv
input/CM_cds.fa
input/CZ_cds.fa
input/CM_proteins_cleaned.faa
input/CZ_proteins_cleaned.faa
input/effector_one_to_one_list.txt
input/CWDE_one_to_one_list.txt
```

The ortholog table must contain columns specified in `config.yaml`, for example:

```text
orthogroup
cm_protein
cz_protein
```

## Usage

Edit `config.yaml` first:

```bash
nano config.yaml
```

Run the full pipeline:

```bash
bash run_full_dnds_pipeline.sh
```

## Outputs

Main output files:

```text
results/id_mapping.tsv
results/ortholog_pairs_clean.tsv
results/alignment_length_stats.csv
results/dnds_enhanced_complete.csv
results/dnds_enhanced_clean.csv
results/accelerated_genes.csv
results/evolutionary_analysis/
```

## Quality control

The default quality filters include:

- minimum effective codon length
- maximum gap percentage
- minimum sequence identity
- removal of internal stop codons
- dS range filtering
- omega outlier filtering

These thresholds can be modified in `config.yaml`.

## Statistical analysis

Evolutionary rates were compared among candidate effectors, CWDEs, and background genes. Wilcoxon rank-sum tests were used for pairwise comparisons. Fisher's exact tests were used to test enrichment of high-omega genes. Multivariable linear regression was used to assess gene category effects while controlling for effective codon length, sequence identity, GC content, and gap percentage.
