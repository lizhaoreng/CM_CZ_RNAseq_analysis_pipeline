#!/usr/bin/env Rscript

# ==============================================================================
# Script: plot_CWDE_expression_heatmap.R
#
# Purpose:
#   1. Read CWDE expression table from an Excel file.
#   2. Automatically identify expression columns.
#   3. Merge multi-copy genes by Orthogroup | CAZy_family, summing expression values.
#   4. Filter CWDE genes with low expression in in vitro and 0 dpi samples but
#      high expression at 3 dpi and/or 6 dpi.
#   5. Export filtered gene tables.
#   6. Draw a fire-style heatmap using ComplexHeatmap.
#
# Expected input:
#   An Excel file containing at least the following columns:
#     - Orthogroup
#     - CAZy_family
#     - CAZy_function
#     - Presence_type
#     - CM_in vitro
#     - CM_0dpi
#     - CM_1dpi
#     - CM_3dpi
#     - CM_6dpi
#     - CZ_in vitro
#     - CZ_0dpi
#     - CZ_1dpi
#     - CZ_3dpi
#     - CZ_6dpi
#
# Optional columns:
#     - CAZy_substrate
#     - CWDE_Quality_Score
#     - High_Level_Category
#
# Usage:
#   Rscript plot_CWDE_expression_heatmap.R
#
# Or modify the parameter section below.
#
# Author: Zhaoreng Li et al.
# ==============================================================================


# ==============================================================================
# 0. Load packages
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
  library(openxlsx)
})


# ==============================================================================
# 1. Parameters
# ==============================================================================

# Input file
input_file <- "CWDE_with_FPKM.xlsx"

# Output directory
output_dir <- "CWDE_heatmap_output"

# Expression column order
desired_expr_order <- c(
  "CM_in vitro", "CM_0dpi", "CM_1dpi", "CM_3dpi", "CM_6dpi",
  "CZ_in vitro", "CZ_0dpi", "CZ_1dpi", "CZ_3dpi", "CZ_6dpi"
)

# Filtering thresholds
early_low_threshold <- 5
late_high_threshold <- 20
fold_change_threshold <- 2

# Maximum number of genes displayed in heatmap
max_genes_for_heatmap <- 100

# Heatmap output size
fig_width <- 10
fig_height <- 10

# Heatmap output resolution
tiff_res <- 600


# ==============================================================================
# 2. Create output directory
# ==============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ==============================================================================
# 3. Helper functions
# ==============================================================================

message_header <- function(x) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(x, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
}


sort_expression_columns <- function(cols) {
  desired_expr_order[desired_expr_order %in% cols]
}


convert_to_matrix_with_rownames <- function(data, rowname_col) {
  row_names <- data[[rowname_col]]
  data_matrix <- data %>%
    dplyr::select(-dplyr::all_of(rowname_col)) %>%
    as.matrix()
  rownames(data_matrix) <- row_names
  return(data_matrix)
}


create_presence_type_colors <- function(presence_types) {
  predefined_colors <- c(
    "Single_Copy" = "#4E79A7",
    "Multi_Copy" = "#E15759",
    "CM_Single_CZ_Multi" = "#59A14F",
    "CM_Multi_CZ_Single" = "#F28E2B",
    "CM_Specific" = "#B07AA1",
    "CZ_Specific" = "#9C755F"
  )
  
  presence_types <- unique(as.character(presence_types))
  result_colors <- predefined_colors[presence_types]
  result_colors[is.na(result_colors)] <- "#BDBDBD"
  
  return(result_colors)
}


create_cazy_function_colors <- function(cazy_functions) {
  predefined_colors <- c(
    "Hemicellulase" = "#FF9999",
    "Cellulase" = "#99FF99",
    "Pectinase" = "#9999FF",
    "Cutinase" = "#FFFF99",
    "Chitinase" = "#FFB6C1"
  )
  
  cazy_functions <- unique(as.character(cazy_functions))
  result_colors <- predefined_colors[cazy_functions]
  
  default_colors <- c(
    "#FFB6C1", "#98FB98", "#87CEEB", "#F0E68C",
    "#DDA0DD", "#B0C4DE", "#FFA07A", "#90EE90"
  )
  
  missing_idx <- which(is.na(result_colors))
  if (length(missing_idx) > 0) {
    result_colors[missing_idx] <- default_colors[
      ((seq_along(missing_idx) - 1) %% length(default_colors)) + 1
    ]
  }
  
  names(result_colors) <- cazy_functions
  return(result_colors)
}


# ==============================================================================
# 4. Read input data
# ==============================================================================

message_header("Reading CWDE expression table")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

cwde_data <- readxl::read_excel(input_file)

cat("Input data dimensions:", paste(dim(cwde_data), collapse = " x "), "\n")
cat("Column names:\n")
print(colnames(cwde_data))

expression_cols <- grep(
  "^(CM|CZ)_(in vitro|0dpi|1dpi|3dpi|6dpi)$",
  colnames(cwde_data),
  value = TRUE
)
expression_cols <- sort_expression_columns(expression_cols)

cat("\nDetected expression columns:\n")
print(expression_cols)

if (length(expression_cols) == 0) {
  stop("No expression columns were detected. Please check column names.")
}


# ==============================================================================
# 5. Heatmap color scale
# ==============================================================================

color_breaks <- c(0, 2, 4, 6, 8)

col_fun_fire <- circlize::colorRamp2(
  color_breaks,
  c("#000080", "#00BFFF", "#FFFFFF", "#FFD700", "#FF4500")
)


# ==============================================================================
# 6. Process CWDE expression data
# ==============================================================================

process_cwde_expression_data <- function(data) {
  cat("Processing CWDE expression data...\n")
  
  required_cols <- c("Orthogroup", "CAZy_family", "CAZy_function", "Presence_type")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  available_expr_cols <- grep(
    "^(CM|CZ)_(in vitro|0dpi|1dpi|3dpi|6dpi)$",
    colnames(data),
    value = TRUE
  )
  available_expr_cols <- sort_expression_columns(available_expr_cols)
  
  if (length(available_expr_cols) == 0) {
    stop("No expression columns found.")
  }
  
  has_substrate <- "CAZy_substrate" %in% colnames(data)
  has_score <- "CWDE_Quality_Score" %in% colnames(data)
  has_category <- "High_Level_Category" %in% colnames(data)
  
  processed_data <- data %>%
    dplyr::filter(!is.na(Orthogroup), Orthogroup != "") %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(available_expr_cols),
        ~ suppressWarnings(as.numeric(as.character(.x)))
      ),
      dplyr::across(
        dplyr::all_of(available_expr_cols),
        ~ ifelse(is.na(.x), 0, .x)
      ),
      Orthogroup_CAZy = paste(Orthogroup, CAZy_family, sep = " | ")
    )
  
  duplicate_check <- processed_data %>%
    dplyr::count(Orthogroup_CAZy) %>%
    dplyr::filter(n > 1)
  
  if (nrow(duplicate_check) > 0) {
    cat("Detected", nrow(duplicate_check),
        "duplicated Orthogroup | CAZy_family combinations.\n")
    cat("Merging multi-copy genes by summing expression values...\n")
    
    grouped_data <- processed_data %>%
      dplyr::group_by(
        Orthogroup_CAZy,
        Orthogroup,
        CAZy_family,
        CAZy_function,
        Presence_type
      ) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(available_expr_cols),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
    
    if (has_substrate) {
      substrate_df <- processed_data %>%
        dplyr::group_by(
          Orthogroup_CAZy,
          Orthogroup,
          CAZy_family,
          CAZy_function,
          Presence_type
        ) %>%
        dplyr::summarise(
          CAZy_substrate = dplyr::first(CAZy_substrate),
          .groups = "drop"
        )
      
      grouped_data <- dplyr::left_join(
        grouped_data,
        substrate_df,
        by = c(
          "Orthogroup_CAZy",
          "Orthogroup",
          "CAZy_family",
          "CAZy_function",
          "Presence_type"
        )
      )
    }
    
    if (has_score) {
      score_df <- processed_data %>%
        dplyr::group_by(
          Orthogroup_CAZy,
          Orthogroup,
          CAZy_family,
          CAZy_function,
          Presence_type
        ) %>%
        dplyr::summarise(
          CWDE_Quality_Score = suppressWarnings(max(CWDE_Quality_Score, na.rm = TRUE)),
          .groups = "drop"
        )
      
      score_df$CWDE_Quality_Score[is.infinite(score_df$CWDE_Quality_Score)] <- NA
      
      grouped_data <- dplyr::left_join(
        grouped_data,
        score_df,
        by = c(
          "Orthogroup_CAZy",
          "Orthogroup",
          "CAZy_family",
          "CAZy_function",
          "Presence_type"
        )
      )
    }
    
    if (has_category) {
      category_df <- processed_data %>%
        dplyr::group_by(
          Orthogroup_CAZy,
          Orthogroup,
          CAZy_family,
          CAZy_function,
          Presence_type
        ) %>%
        dplyr::summarise(
          High_Level_Category = dplyr::first(High_Level_Category),
          .groups = "drop"
        )
      
      grouped_data <- dplyr::left_join(
        grouped_data,
        category_df,
        by = c(
          "Orthogroup_CAZy",
          "Orthogroup",
          "CAZy_family",
          "CAZy_function",
          "Presence_type"
        )
      )
    }
    
    processed_data <- grouped_data
  }
  
  cat("Processed data dimensions:",
      paste(dim(processed_data), collapse = " x "), "\n")
  
  return(processed_data)
}


# ==============================================================================
# 7. Filter CWDE genes
# ==============================================================================

filter_low_early_high_late_genes <- function(data,
                                             early_low_threshold = 5,
                                             late_high_threshold = 20,
                                             fold_change_threshold = 2) {
  cat("Filtering CWDE genes with low early expression and high late expression...\n")
  
  required_cols <- c(
    "CM_in vitro", "CM_0dpi", "CM_3dpi", "CM_6dpi",
    "CZ_in vitro", "CZ_0dpi", "CZ_3dpi", "CZ_6dpi"
  )
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing expression columns required for filtering: ",
         paste(missing_cols, collapse = ", "))
  }
  
  filtered_data <- data %>%
    dplyr::mutate(
      CM_Early_Mean = (`CM_in vitro` + CM_0dpi) / 2,
      CM_Late_Mean  = (CM_3dpi + CM_6dpi) / 2,
      CZ_Early_Mean = (`CZ_in vitro` + CZ_0dpi) / 2,
      CZ_Late_Mean  = (CZ_3dpi + CZ_6dpi) / 2,
      
      CM_Late_vs_Early_FC = (CM_Late_Mean + 1) / (CM_Early_Mean + 1),
      CZ_Late_vs_Early_FC = (CZ_Late_Mean + 1) / (CZ_Early_Mean + 1),
      
      CM_Pass = (CM_Early_Mean <= early_low_threshold) &
        (CM_Late_Mean >= late_high_threshold) &
        (CM_Late_vs_Early_FC >= fold_change_threshold),
      
      CZ_Pass = (CZ_Early_Mean <= early_low_threshold) &
        (CZ_Late_Mean >= late_high_threshold) &
        (CZ_Late_vs_Early_FC >= fold_change_threshold),
      
      Pass_Either = CM_Pass | CZ_Pass,
      Pass_Both = CM_Pass & CZ_Pass
    ) %>%
    dplyr::filter(Pass_Either)
  
  cat("Number of filtered CWDE genes:", nrow(filtered_data), "\n")
  
  return(filtered_data)
}


# ==============================================================================
# 8. Draw fire-style CWDE heatmap
# ==============================================================================

create_cwde_heatmap_fire <- function(data,
                                     max_genes = 100,
                                     cluster_rows = TRUE,
                                     output_dir = ".",
                                     fig_width = 10,
                                     fig_height = 10,
                                     tiff_res = 600) {
  cat("Drawing fire-style CWDE heatmap...\n")
  
  available_expr_cols <- grep(
    "^(CM|CZ)_(in vitro|0dpi|1dpi|3dpi|6dpi)$",
    colnames(data),
    value = TRUE
  )
  available_expr_cols <- sort_expression_columns(available_expr_cols)
  
  if (length(available_expr_cols) == 0) {
    stop("No expression columns found for heatmap.")
  }
  
  expr_data <- data %>%
    dplyr::filter(!is.na(Orthogroup_CAZy)) %>%
    dplyr::rowwise() %>%
    dplyr::filter(
      sum(dplyr::c_across(dplyr::all_of(available_expr_cols)), na.rm = TRUE) > 0
    ) %>%
    dplyr::ungroup()
  
  if (nrow(expr_data) == 0) {
    stop("No valid rows available for heatmap.")
  }
  
  if (nrow(expr_data) > max_genes) {
    expr_data <- expr_data %>%
      dplyr::mutate(
        Heatmap_Score = (
          (CM_3dpi + CM_6dpi + CZ_3dpi + CZ_6dpi) -
            (`CM_in vitro` + CM_0dpi + `CZ_in vitro` + CZ_0dpi)
        )
      ) %>%
      dplyr::arrange(desc(Heatmap_Score)) %>%
      dplyr::slice_head(n = max_genes)
    
    cat("Heatmap retained top", max_genes,
        "genes with the strongest late-induction pattern.\n")
  }
  
  expr_matrix_data <- expr_data %>%
    dplyr::select(Orthogroup_CAZy, dplyr::all_of(available_expr_cols))
  
  expr_matrix <- convert_to_matrix_with_rownames(expr_matrix_data, "Orthogroup_CAZy")
  expr_matrix[is.na(expr_matrix)] <- 0
  expr_matrix[is.infinite(expr_matrix)] <- 0
  
  expr_matrix_log <- log2(expr_matrix + 1)
  
  row_annotation_df <- expr_data %>%
    dplyr::select(Orthogroup_CAZy, Presence_type, CAZy_function) %>%
    as.data.frame()
  
  rownames(row_annotation_df) <- row_annotation_df$Orthogroup_CAZy
  row_annotation_df <- row_annotation_df[rownames(expr_matrix_log), , drop = FALSE]
  
  annotation_colors <- list(
    Presence_type = create_presence_type_colors(row_annotation_df$Presence_type),
    CAZy_function = create_cazy_function_colors(row_annotation_df$CAZy_function)
  )
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Presence_type = row_annotation_df$Presence_type,
    CAZy_function = row_annotation_df$CAZy_function,
    col = annotation_colors,
    annotation_name_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold.italic"
    ),
    annotation_name_rot = 45,
    simple_anno_size = grid::unit(0.5, "cm"),
    annotation_name_side = "bottom",
    gap = grid::unit(2, "mm")
  )
  
  col_info <- data.frame(
    Sample = colnames(expr_matrix_log),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Strain = ifelse(grepl("^CM_", Sample), "CM", "CZ"),
      Time_Point = sub("^(CM|CZ)_(in vitro|0dpi|1dpi|3dpi|6dpi)$", "\\2", Sample)
    )
  
  rownames(col_info) <- col_info$Sample
  
  time_point_colors <- c(
    "in vitro" = "#D9EAD3",
    "0dpi"  = "#CFE2F3",
    "1dpi"  = "#FFF2CC",
    "3dpi"  = "#FFE599",
    "6dpi"  = "#FFD966"
  )
  
  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    Strain = col_info$Strain,
    Time_Point = col_info$Time_Point,
    col = list(
      Strain = c("CM" = "#FF6B6B", "CZ" = "#4ECDC4"),
      Time_Point = time_point_colors
    ),
    annotation_name_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
    annotation_name_rot = 0,
    simple_anno_size = grid::unit(0.4, "cm"),
    gap = grid::unit(2, "mm")
  )
  
  ht <- ComplexHeatmap::Heatmap(
    expr_matrix_log,
    name = "Log2(FPKM+1)",
    col = col_fun_fire,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    cluster_rows = cluster_rows,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 45,
    left_annotation = row_ha,
    top_annotation = col_ha,
    column_title = NULL,
    width = grid::unit(12, "cm"),
    height = grid::unit(max(8, nrow(expr_matrix_log) * 0.22), "cm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.5)
  )
  
  pdf_file <- file.path(
    output_dir,
    "CWDE_heatmap_fire_filtered_low_early_high_late.pdf"
  )
  
  tiff_file <- file.path(
    output_dir,
    "CWDE_heatmap_fire_filtered_low_early_high_late.tiff"
  )
  
  pdf(pdf_file, width = fig_width, height = fig_height)
  ComplexHeatmap::draw(ht)
  dev.off()
  
  tiff(
    filename = tiff_file,
    width = fig_width,
    height = fig_height,
    units = "in",
    res = tiff_res,
    compression = "lzw"
  )
  ComplexHeatmap::draw(ht)
  dev.off()
  
  cat("Heatmap files saved:\n")
  cat(" -", pdf_file, "\n")
  cat(" -", tiff_file, "\n")
  cat("Number of genes shown in heatmap:", nrow(expr_matrix_log), "\n")
  
  return(ht)
}


# ==============================================================================
# 9. Run analysis
# ==============================================================================

message_header("Data preprocessing")

processed_cwde_data <- process_cwde_expression_data(cwde_data)

cat("\nPresence_type distribution:\n")
print(table(processed_cwde_data$Presence_type))

cat("\nCAZy_function distribution:\n")
print(table(processed_cwde_data$CAZy_function))


message_header("Filtering target CWDE genes")

filtered_cwde_data <- filter_low_early_high_late_genes(
  processed_cwde_data,
  early_low_threshold = early_low_threshold,
  late_high_threshold = late_high_threshold,
  fold_change_threshold = fold_change_threshold
)

if (nrow(filtered_cwde_data) == 0) {
  stop("No genes passed the filtering criteria. Please use less stringent thresholds.")
}


message_header("Exporting filtered gene tables")

filtered_output <- filtered_cwde_data %>%
  dplyr::arrange(
    desc(Pass_Both),
    desc(CM_Late_vs_Early_FC + CZ_Late_vs_Early_FC)
  )

csv_file <- file.path(
  output_dir,
  "CWDE_filtered_low_early_high_late.csv"
)

xlsx_file <- file.path(
  output_dir,
  "CWDE_filtered_low_early_high_late.xlsx"
)

write.csv(
  filtered_output,
  csv_file,
  row.names = FALSE,
  quote = FALSE
)

openxlsx::write.xlsx(
  filtered_output,
  xlsx_file,
  rowNames = FALSE
)

cat("Filtered gene tables saved:\n")
cat(" -", csv_file, "\n")
cat(" -", xlsx_file, "\n")


message_header("Drawing heatmap")

cwde_heatmap_fire <- create_cwde_heatmap_fire(
  filtered_cwde_data,
  max_genes = max_genes_for_heatmap,
  cluster_rows = TRUE,
  output_dir = output_dir,
  fig_width = fig_width,
  fig_height = fig_height,
  tiff_res = tiff_res
)


# ==============================================================================
# 10. Save session information
# ==============================================================================

session_file <- file.path(output_dir, "sessionInfo_CWDE_heatmap.txt")

sink(session_file)
cat("Session information\n")
cat("===================\n\n")
print(sessionInfo())
sink()

cat("Session information saved:\n")
cat(" -", session_file, "\n")


# ==============================================================================
# 11. Analysis summary
# ==============================================================================

message_header("CWDE filtering and heatmap analysis completed")

cat("Filtering criteria:\n")
cat(sprintf(" - Early low-expression threshold: mean(in vitro, 0 dpi) <= %s\n",
            early_low_threshold))
cat(sprintf(" - Late high-expression threshold: mean(3 dpi, 6 dpi) >= %s\n",
            late_high_threshold))
cat(sprintf(" - Late/Early fold-change threshold: >= %s\n",
            fold_change_threshold))

cat("\nData summary:\n")
cat(sprintf(" - Total processed CWDE entries: %d\n", nrow(processed_cwde_data)))
cat(sprintf(" - Filtered CWDE entries: %d\n", nrow(filtered_cwde_data)))
cat(sprintf(" - CM-passed entries: %d\n", sum(filtered_cwde_data$CM_Pass, na.rm = TRUE)))
cat(sprintf(" - CZ-passed entries: %d\n", sum(filtered_cwde_data$CZ_Pass, na.rm = TRUE)))
cat(sprintf(" - Entries passing both CM and CZ criteria: %d\n",
            sum(filtered_cwde_data$Pass_Both, na.rm = TRUE)))

cat("\nColumn order used in heatmap:\n")
for (x in desired_expr_order) {
  cat(sprintf(" - %s\n", x))
}

cat("\nOutput files:\n")
cat(" -", csv_file, "\n")
cat(" -", xlsx_file, "\n")
cat(" -", file.path(output_dir, "CWDE_heatmap_fire_filtered_low_early_high_late.pdf"), "\n")
cat(" -", file.path(output_dir, "CWDE_heatmap_fire_filtered_low_early_high_late.tiff"), "\n")
cat(" -", session_file, "\n")

cat("\nHeatmap notes:\n")
cat(" - Rows: filtered CWDE genes or Orthogroup | CAZy_family entries\n")
cat(" - Right row labels: Orthogroup_CAZy\n")
cat(" - Left row annotations: Presence_type and CAZy_function\n")
cat(" - Row annotation titles: 45 degrees, bold italic\n")
cat(" - Top column annotations: Strain and Time_Point, not italicized\n")
cat(" - Heatmap legend title: Log2(FPKM+1)\n")
cat(" - Heatmap title: removed\n")
cat(" - Color scheme: fire-style color scale\n")

cat("\nAnalysis completed successfully.\n")
