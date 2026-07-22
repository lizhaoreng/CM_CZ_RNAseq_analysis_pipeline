#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_run_deseq2_pairwise_comparisons.R
#
# Purpose:
#   Run DESeq2 pairwise comparisons for CM and CZ pathogen-derived count matrices.
#   The script compares post-inoculation samples with 0 dpi and in vitro samples.
#
# Input files:
#   - CM_counts.csv
#   - CZ_counts.csv
#   - CM_FPKM.csv
#   - CZ_FPKM.csv
#
# Output files:
#   - results/deseq2/CM_DESeq2_result.csv
#   - results/deseq2/CZ_DESeq2_result.csv
#
# Required packages:
#   DESeq2, dplyr, readr
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
})

# -----------------------------
# Parameters
# -----------------------------

input_dir <- "data"
output_dir <- "results/deseq2"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_config <- list(
  CM = list(
    counts_file = file.path(input_dir, "CM_counts.csv"),
    fpkm_file   = file.path(input_dir, "CM_FPKM.csv"),
    prefix      = "CM"
  ),
  CZ = list(
    counts_file = file.path(input_dir, "CZ_counts.csv"),
    fpkm_file   = file.path(input_dir, "CZ_FPKM.csv"),
    prefix      = "CZ"
  )
)

comparisons <- list(
  "1dpi_vs_0dpi"     = c("condition", "1dpi", "0dpi"),
  "3dpi_vs_0dpi"     = c("condition", "3dpi", "0dpi"),
  "6dpi_vs_0dpi"     = c("condition", "6dpi", "0dpi"),
  "1dpi_vs_invitro"  = c("condition", "1dpi", "in vitro"),
  "3dpi_vs_invitro"  = c("condition", "3dpi", "in vitro"),
  "6dpi_vs_invitro"  = c("condition", "6dpi", "in vitro")
)

# -----------------------------
# Helper functions
# -----------------------------

check_file <- function(file) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }
}

infer_condition <- function(sample_names, prefix) {
  condition <- sub(paste0("^", prefix, "_"), "", sample_names)
  condition <- sub("_[0-9]+$", "", condition)
  condition
}

run_deseq2_comparisons <- function(counts_file, fpkm_file, prefix, comparisons) {
  check_file(counts_file)
  check_file(fpkm_file)
  
  message("\nRunning DESeq2 for ", prefix)
  
  counts_raw <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
  fpkm_raw   <- read.csv(fpkm_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  if (!"Geneid" %in% colnames(counts_raw)) {
    stop("Column 'Geneid' was not found in ", counts_file)
  }
  
  if (!"Geneid" %in% colnames(fpkm_raw)) {
    stop("Column 'Geneid' was not found in ", fpkm_file)
  }
  
  rownames(counts_raw) <- counts_raw$Geneid
  counts_mat <- as.matrix(counts_raw[, setdiff(colnames(counts_raw), "Geneid"), drop = FALSE])
  storage.mode(counts_mat) <- "integer"
  
  rownames(fpkm_raw) <- fpkm_raw$Geneid
  
  sample_names <- colnames(counts_mat)
  condition <- infer_condition(sample_names, prefix)
  
  col_data <- data.frame(
    row.names = sample_names,
    condition = factor(condition)
  )
  
  message("Sample grouping:")
  print(table(col_data$condition))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = col_data,
    design = ~ condition
  )
  
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds)
  
  result_list <- lapply(names(comparisons), function(comp_name) {
    contrast <- comparisons[[comp_name]]
    
    if (!all(contrast[2:3] %in% levels(col_data$condition))) {
      warning("Skipping comparison because condition level is missing: ", comp_name)
      return(NULL)
    }
    
    res <- results(dds, contrast = contrast, independentFiltering = TRUE)
    res_df <- as.data.frame(res)
    
    res_df <- res_df[, c("log2FoldChange", "padj"), drop = FALSE]
    res_df$Geneid <- rownames(res_df)
    rownames(res_df) <- NULL
    
    colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- paste0("log2FC_", comp_name)
    colnames(res_df)[colnames(res_df) == "padj"] <- paste0("FDR_", comp_name)
    
    res_df[, c("Geneid", paste0("log2FC_", comp_name), paste0("FDR_", comp_name))]
  })
  
  result_list <- result_list[!sapply(result_list, is.null)]
  
  deg_wide <- Reduce(function(a, b) merge(a, b, by = "Geneid", all = TRUE), result_list)
  
  counts_df <- counts_raw
  rownames(counts_df) <- NULL
  cnt_cols <- setdiff(colnames(counts_df), "Geneid")
  colnames(counts_df)[colnames(counts_df) %in% cnt_cols] <- paste0("counts_", cnt_cols)
  
  fpkm_df <- fpkm_raw
  rownames(fpkm_df) <- NULL
  fpkm_cols <- setdiff(colnames(fpkm_df), "Geneid")
  colnames(fpkm_df)[colnames(fpkm_df) %in% fpkm_cols] <- paste0("FPKM_", fpkm_cols)
  
  final_df <- merge(counts_df, fpkm_df, by = "Geneid", all.x = TRUE)
  final_df <- merge(final_df, deg_wide, by = "Geneid", all.x = TRUE)
  
  final_df
}

# -----------------------------
# Run analysis
# -----------------------------

for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  
  result <- run_deseq2_comparisons(
    counts_file = cfg$counts_file,
    fpkm_file = cfg$fpkm_file,
    prefix = cfg$prefix,
    comparisons = comparisons
  )
  
  out_file <- file.path(output_dir, paste0(sp, "_DESeq2_result.csv"))
  write.csv(result, out_file, row.names = FALSE, quote = FALSE)
  
  message("Saved: ", out_file)
  message(sp, " result dimensions: ", nrow(result), " rows x ", ncol(result), " columns")
}

session_file <- file.path(output_dir, "sessionInfo_DESeq2.txt")
sink(session_file)
print(sessionInfo())
sink()

message("\nDESeq2 analysis completed.")
