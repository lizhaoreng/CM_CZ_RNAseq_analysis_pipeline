# CM-CZ Comparative  Transcriptomics Pipelines

This repository contains scripts and workflow files used for comparative genomics, pathogen-focused RNA-seq analysis, GO enrichment,  secretome/effector prediction, marker phylogeny, and functional annotation of representative *Cercospora zeae-maydis* and *Cercospora zeina* isolates.

## Directory structure

- `rnaseq_pipeline/`: RNA-seq read processing, mapping, quantification, and differential expression analysis.
- `DeSeq2_GO_pipeline/`: DESeq2 differential expression and GO enrichment analysis.
- `dnds_pipeline/`: pairwise ortholog and dN/dS analysis.
- `functional_annotation_pipeline/`: CAZyme, PHI-base, KEGG, eggNOG, and InterPro annotation scripts.
- `secretome_effector_pipeline/`: secretome and candidate effector prediction.
- `marker_phylogeny_pipeline/`: multilocus phylogenetic analysis using ITS, TEF1, ACT, CAL, and HIS3.
- `heatmap/`: scripts for expression heatmap visualization.

## Notes

Large sequencing files and intermediate alignment files are not included in this repository.
