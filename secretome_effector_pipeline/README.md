# Secretome and effector prediction pipeline

This repository contains scripts used to predict classical secreted proteins and candidate effectors from fungal protein sequences.

## Overview

The pipeline performs:

1. Protein FASTA preparation and format checking.
2. Protein sequence cleaning.
3. Signal peptide prediction using SignalP6.
4. Transmembrane helix filtering using TMHMM.
5. Candidate effector prediction using EffectorP.
6. Subcellular localization prediction of candidate effectors using DeepLoc2.

## Input files

Protein FASTA files should be placed in:

```text
ref/
  CM_proteins_cleaned.faa
  CZ_proteins_cleaned.faa
```

Input file paths can be modified in `config.sh`.

## Usage

Edit the configuration file:

```bash
nano config.sh
```

Run the full pipeline:

```bash
bash run_secretome_effector_pipeline.sh
```

Or run each step separately:

```bash
bash scripts/01_prepare_secretome_input.sh config.sh
bash scripts/02_signalp_tmhmm_filter.sh config.sh
bash scripts/03_effectorp_prediction.sh config.sh
bash scripts/04_deeploc2_effector_localization.sh config.sh
```

## Required software

- Python 3
- Biopython
- SignalP6
- TMHMM
- EffectorP
- DeepLoc2

## Outputs

Main outputs are saved under:

```text
results/secretome/
```

Important files include:

```text
results/secretome/02_signalp_tmhmm_summary.tsv
results/secretome/03_effectorp_summary.tsv
results/secretome/CM/CM_signalp_proteins.faa
results/secretome/CZ/CZ_signalp_proteins.faa
results/secretome/CM/CM_secreted_proteins.faa
results/secretome/CZ/CZ_secreted_proteins.faa
results/secretome/CM/effectorp/CM_effectors.faa
results/secretome/CZ/effectorp/CZ_effectors.faa
results/secretome/CM/deeploc2_effector/
results/secretome/CZ/deeploc2_effector/
```

## Definition of classical secreted proteins

Proteins were considered classical secreted protein candidates if they:

1. Had a predicted signal peptide by SignalP6.
2. Had no predicted transmembrane helix according to TMHMM.

Candidate effectors were then predicted from the classical secreted protein set using EffectorP.
