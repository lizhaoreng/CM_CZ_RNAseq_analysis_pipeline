# Differential expression and GO enrichment analysis

This repository contains R scripts used for pathogen-focused RNA-seq differential expression analysis and GO enrichment analysis.

## Scripts

```text
scripts/
  01_run_deseq2_pairwise_comparisons.R
  02_go_enrichment_infection_up_genes.R
  03_go_enrichment_specific_expressed_genes.R
```

## Required input files

```text
data/
  CM_counts.csv
  CZ_counts.csv
  CM_FPKM.csv
  CZ_FPKM.csv
  background_GO.xlsx
  CM_infection_UP.csv
  CZ_infection_UP.csv
  CM_specific.csv
  CZ_specific.csv
```

## DESeq2 analysis

The script `01_run_deseq2_pairwise_comparisons.R` performs pairwise comparisons:

- 1 dpi vs 0 dpi
- 3 dpi vs 0 dpi
- 6 dpi vs 0 dpi
- 1 dpi vs in vitro
- 3 dpi vs in vitro
- 6 dpi vs in vitro

Run:

```bash
Rscript scripts/01_run_deseq2_pairwise_comparisons.R
```

Outputs:

```text
results/deseq2/CM_DESeq2_result.csv
results/deseq2/CZ_DESeq2_result.csv
```

## GO enrichment for infection-upregulated genes

Run:

```bash
Rscript scripts/02_go_enrichment_infection_up_genes.R
```

Outputs:

```text
results/go_infection_up/
```

## GO enrichment for specific expressed genes

Run:

```bash
Rscript scripts/03_go_enrichment_specific_expressed_genes.R
```

Outputs:

```text
results/go_specific_expressed/
```

## Required R packages

- DESeq2
- dplyr
- tidyr
- readr
- readxl
- stringr
- ggplot2
- GO.db
- tidytext

## Notes

GO annotations from eggNOG and InterPro were combined for enrichment analysis. Fisher's exact test was used to test GO term enrichment, and Benjamini-Hochberg correction was applied to control the false discovery rate. Broad and non-informative GO terms were removed before visualization.
