#!/usr/bin/env Rscript

# ==============================================================================
# Script: plot_effector_expression_heatmap.R
#
# Purpose:
#   1. Read an effector orthogroup expression and annotation table.
#   2. Support both mean-expression columns and replicate-level FPKM columns.
#   3. Automatically calculate mean FPKM values from three biological replicates.
#   4. Filter candidate effector orthogroups showing low expression in early
#      conditions and high expression at later post-inoculation stages.
#   5. Force-retain selected representative effector or pathogenicity-related
#      orthogroups specified in `special_genes`.
#   6. Draw a fire-style ComplexHeatmap with Presence_type row annotation.
#   7. Export filtered tables, heatmap figures, analysis report, and session info.
#
# Expected input:
#   An Excel table containing an orthogroup column and expression columns.
#
# Required identifier column:
#   - orthogroup
#     or one of similar names containing "Orthogroup", "orthogroup", or "Ortho"
#
# Supported mean-expression columns:
#   - CM_in vitro
#   - CM_0dpi
#   - CM_1dpi
#   - CM_3dpi
#   - CM_6dpi
#   - CZ_in vitro
#   - CZ_0dpi
#   - CZ_1dpi
#   - CZ_3dpi
#   - CZ_6dpi
#
# Supported replicate-level columns:
#   - CM_in vitro_1, CM_in vitro_2, CM_in vitro_3
#   - CM_0dpi_1,     CM_0dpi_2,     CM_0dpi_3
#   - ...
#   - CZ_6dpi_1,     CZ_6dpi_2,     CZ_6dpi_3
#
# Also supports old-style prefixes:
#   - CM_expr_CM_in vitro_1 -> CM_in vitro_1
#   - CZ_expr_CZ_3dpi_2    -> CZ_3dpi_2
#   - CM_expr_in vitro_1   -> CM_in vitro_1
#   - CZ_expr_3dpi_2       -> CZ_3dpi_2
#
# Optional annotation columns:
#   - Presence_type
#   - Effector_status
#   - Relation_Type
#   - Copy_Type
#   - copy_type
#   - CM_Gene.x and CZ_Gene.x
#
# Usage:
#   Rscript plot_effector_expression_heatmap.R
#
# Author: Zhaoreng Li et al.
# ==============================================================================


# ==============================================================================
# 0. Load packages
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
  library(openxlsx)
  library(stringr)
})


# ==============================================================================
# 1. Parameters
# ==============================================================================

# Candidate input files. The first existing file will be used.
candidate_files <- c(
  "Effector_OG_expression_annotation_all.xlsx",
  "All_gene_anno_reclassified_annotations.PHI_target_fixed.simple.xlsx",
  "All_gene_anno.xlsx"
)

# Output directory
output_dir <- "effector_heatmap_output"

# Expression column order
desired_expr_order <- c(
  "CM_in vitro", "CM_0dpi", "CM_1dpi", "CM_3dpi", "CM_6dpi",
  "CZ_in vitro", "CZ_0dpi", "CZ_1dpi", "CZ_3dpi", "CZ_6dpi"
)

# Filtering thresholds
early_low_threshold <- 5
late_high_threshold <- 20
fold_change_threshold <- 2

# Maximum number of genes shown in heatmap
max_genes_for_heatmap <- 150

# Heatmap output size
fig_width <- 10
fig_height <- 8

# TIFF resolution
tiff_res <- 600

# Representative effector/pathogenicity-related genes to retain in the heatmap
# even if they do not meet the filtering criteria.
special_genes <- c(
  "OG0002164" = "Ecp6-like",
  "OG0000299" = "Ecp2-like",
  "OG0000492" = "Ecp20-2",
  "OG0003305" = "NPP1",
  "OG0008208" = "Nis1-like",
  "OG0001101" = "PsXEG1",
  "OG0006941" = "SsCut1-like",
  "OG0000107" = "BGL2-like",
  "OG0002278" = "Fpr2",
  "OG0006905" = "Flc/Pkd2-like"
)


# ==============================================================================
# 2. Create output directory and select input file
# ==============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

input_file <- candidate_files[file.exists(candidate_files)][1]

if (is.na(input_file) || length(input_file) == 0) {
  stop(
    "No input file was found. Please provide at least one of the following files:\n",
    paste(candidate_files, collapse = "\n")
  )
}


# ==============================================================================
# 3. Helper functions
# ==============================================================================

message_header <- function(x) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(x, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
}


clean_id_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- stringr::str_trim(x)
  return(x)
}


safe_numeric <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "nan", "NULL", "null", "None", "-")] <- NA
  suppressWarnings(as.numeric(x))
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


make_display_label <- function(orthogroup_ids, special_gene_map) {
  sapply(
    orthogroup_ids,
    function(x) {
      if (x %in% names(special_gene_map)) {
        paste0(x, " | ", special_gene_map[[x]])
      } else {
        x
      }
    },
    USE.NAMES = FALSE
  )
}


create_presence_type_colors <- function(type_values) {
  predefined_colors <- c(
    "Single_Copy" = "#4E79A7",
    "Multi_Copy" = "#E15759",
    "CM_Single_CZ_Multi" = "#59A14F",
    "CM_Multi_CZ_Single" = "#F28E2B",
    "CM_Specific" = "#B07AA1",
    "CZ_Specific" = "#9C755F",
    "Both_effectors" = "#1F77B4",
    "CM_effector_only" = "#D62728",
    "CZ_effector_only" = "#2CA02C",
    "Both_Single" = "#4E79A7",
    "Unknown" = "#BDBDBD"
  )
  
  type_values <- unique(as.character(type_values))
  result_colors <- predefined_colors[type_values]
  result_colors[is.na(result_colors)] <- "#BDBDBD"
  names(result_colors) <- type_values
  
  return(result_colors)
}


# ==============================================================================
# 4. Heatmap color scale
# ==============================================================================

color_breaks <- c(0, 2, 4, 6, 8)

col_fun_fire <- circlize::colorRamp2(
  color_breaks,
  c("#000080", "#00BFFF", "#FFFFFF", "#FFD700", "#FF4500")
)


# ==============================================================================
# 5. Read input data
# ==============================================================================

message_header("Reading effector expression table")

cat("Input file:", input_file, "\n")

effector_data <- readxl::read_excel(input_file)

cat("Input data dimensions:", paste(dim(effector_data), collapse = " x "), "\n")
cat("Column names:\n")
print(colnames(effector_data))

cat("\nPreview of first five rows:\n")
print(utils::head(effector_data, 5))


# ==============================================================================
# 6. Standardize expression column names
# ==============================================================================

standardize_expression_column_names <- function(data) {
  cat("\nStandardizing expression column names...\n")
  
  original_names <- colnames(data)
  new_names <- original_names
  
  # Old integrated-table formats:
  # CM_expr_CM_in vitro_1 -> CM_in vitro_1
  # CZ_expr_CZ_3dpi_2    -> CZ_3dpi_2
  new_names <- gsub("^CM_expr_CM_", "CM_", new_names)
  new_names <- gsub("^CZ_expr_CZ_", "CZ_", new_names)
  
  # Additional possible formats:
  # CM_expr_in vitro_1 -> CM_in vitro_1
  # CZ_expr_3dpi_2     -> CZ_3dpi_2
  new_names <- gsub("^CM_expr_", "CM_", new_names)
  new_names <- gsub("^CZ_expr_", "CZ_", new_names)
  
  colnames(data) <- new_names
  
  changed <- data.frame(
    old = original_names,
    new = new_names,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(old != new)
  
  if (nrow(changed) > 0) {
    cat("The following column names were standardized:\n")
    print(changed)
  } else {
    cat("No expression column prefixes required standardization.\n")
  }
  
  return(data)
}


# ==============================================================================
# 7. Calculate mean expression columns from replicate columns
# ==============================================================================

add_mean_expression_columns <- function(data) {
  cat("\nCalculating mean expression columns from replicate-level FPKM columns...\n")
  
  time_points <- c("in vitro", "0dpi", "1dpi", "3dpi", "6dpi")
  strains <- c("CM", "CZ")
  
  for (strain in strains) {
    for (tp in time_points) {
      mean_col <- paste0(strain, "_", tp)
      rep_cols <- paste0(strain, "_", tp, "_", 1:3)
      rep_cols_present <- rep_cols[rep_cols %in% colnames(data)]
      
      # If mean column already exists, keep it after numeric conversion.
      if (mean_col %in% colnames(data)) {
        data[[mean_col]] <- safe_numeric(data[[mean_col]])
        data[[mean_col]][is.na(data[[mean_col]])] <- 0
        cat("Existing mean column retained:", mean_col, "\n")
        next
      }
      
      # If replicate columns are available, calculate row means.
      if (length(rep_cols_present) > 0) {
        for (rc in rep_cols_present) {
          data[[rc]] <- safe_numeric(data[[rc]])
        }
        
        data[[mean_col]] <- rowMeans(
          data[, rep_cols_present, drop = FALSE],
          na.rm = TRUE
        )
        
        data[[mean_col]][is.nan(data[[mean_col]])] <- 0
        data[[mean_col]][is.na(data[[mean_col]])] <- 0
        
        cat(
          "Mean column calculated:", mean_col,
          "from", paste(rep_cols_present, collapse = ", "), "\n"
        )
        
      } else {
        cat("Warning: replicate columns for", mean_col,
            "were not found. This column may be required later.\n")
      }
    }
  }
  
  return(data)
}


# ==============================================================================
# 8. Process effector expression data
# ==============================================================================

process_effector_expression_data <- function(data) {
  cat("\nProcessing effector expression data...\n")
  
  data <- standardize_expression_column_names(data)
  
  # Standardize orthogroup column.
  if (!"orthogroup" %in% colnames(data)) {
    possible_orthogroup_cols <- grep(
      "^orthogroup$|Orthogroup|orthogroup|Ortho|ortho",
      colnames(data),
      value = TRUE,
      ignore.case = TRUE
    )
    
    if (length(possible_orthogroup_cols) > 0) {
      data <- data %>%
        dplyr::rename(orthogroup = !!possible_orthogroup_cols[1])
      
      cat("Column", possible_orthogroup_cols[1],
          "renamed to orthogroup.\n")
    } else {
      stop("No orthogroup column was found.")
    }
  }
  
  data <- add_mean_expression_columns(data)
  
  expr_cols <- grep(
    "^(CM|CZ)_(in vitro|0dpi|1dpi|3dpi|6dpi)$",
    colnames(data),
    value = TRUE
  )
  expr_cols <- sort_expression_columns(expr_cols)
  
  cat("\nDetected mean-expression columns:\n")
  print(expr_cols)
  
  if (length(expr_cols) == 0) {
    stop(
      "No expression columns found. Please provide either mean-expression columns ",
      "or replicate-level columns for CM/CZ in vitro, 0 dpi, 1 dpi, 3 dpi, and 6 dpi."
    )
  }
  
  missing_expr <- setdiff(desired_expr_order, expr_cols)
  if (length(missing_expr) > 0) {
    cat("Warning: the following expected expression columns are missing:\n")
    print(missing_expr)
  }
  
  processed_data <- data %>%
    dplyr::mutate(orthogroup = clean_id_text(orthogroup)) %>%
    dplyr::filter(!is.na(orthogroup), orthogroup != "") %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(expr_cols), ~ safe_numeric(.x)),
      dplyr::across(dplyr::all_of(expr_cols), ~ ifelse(is.na(.x), 0, .x))
    )
  
  # Determine Presence_type.
  if ("Presence_type" %in% colnames(processed_data)) {
    cat("Presence_type column detected and used as row annotation.\n")
    
  } else if ("Effector_status" %in% colnames(processed_data)) {
    processed_data <- processed_data %>%
      dplyr::mutate(Presence_type = Effector_status)
    cat("Presence_type was generated from Effector_status.\n")
    
  } else if ("Relation_Type" %in% colnames(processed_data)) {
    processed_data <- processed_data %>%
      dplyr::mutate(Presence_type = Relation_Type)
    cat("Presence_type was generated from Relation_Type.\n")
    
  } else if ("Copy_Type" %in% colnames(processed_data)) {
    processed_data <- processed_data %>%
      dplyr::mutate(Presence_type = Copy_Type)
    cat("Presence_type was generated from Copy_Type.\n")
    
  } else if ("copy_type" %in% colnames(processed_data)) {
    processed_data <- processed_data %>%
      dplyr::mutate(Presence_type = copy_type)
    cat("Presence_type was generated from copy_type.\n")
    
  } else if ("CM_Gene.x" %in% colnames(processed_data) &&
             "CZ_Gene.x" %in% colnames(processed_data)) {
    processed_data <- processed_data %>%
      dplyr::mutate(
        CM_Gene.x = clean_id_text(CM_Gene.x),
        CZ_Gene.x = clean_id_text(CZ_Gene.x),
        Presence_type = dplyr::case_when(
          CM_Gene.x != "" & CZ_Gene.x != "" ~ "Both_Single",
          CM_Gene.x != "" & CZ_Gene.x == "" ~ "CM_Specific",
          CM_Gene.x == "" & CZ_Gene.x != "" ~ "CZ_Specific",
          TRUE ~ "Unknown"
        )
      )
    cat("Presence_type was generated from CM_Gene.x and CZ_Gene.x.\n")
    
  } else {
    processed_data <- processed_data %>%
      dplyr::mutate(Presence_type = "Unknown")
    cat("No Presence_type-related column was found. Presence_type set to Unknown.\n")
  }
  
  processed_data <- processed_data %>%
    dplyr::mutate(Presence_type = clean_id_text(Presence_type))
  
  return(processed_data)
}


# ==============================================================================
# 9. Filter low-early and high-late expressed effectors
# ==============================================================================

filter_low_early_high_late_effectors <- function(data,
                                                 early_low_threshold = 5,
                                                 late_high_threshold = 20,
                                                 fold_change_threshold = 2,
                                                 special_gene_ids = NULL) {
  cat("\nFiltering effectors with low early expression and high late expression...\n")
  
  required_cols <- c(
    "CM_in vitro", "CM_0dpi", "CM_3dpi", "CM_6dpi",
    "CZ_in vitro", "CZ_0dpi", "CZ_3dpi", "CZ_6dpi"
  )
  
  missing_required_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_required_cols) > 0) {
    cat("Available columns:\n")
    print(colnames(data))
    stop(
      "Missing expression columns required for filtering: ",
      paste(missing_required_cols, collapse = ", ")
    )
  }
  
  if (is.null(special_gene_ids)) {
    special_gene_ids <- character(0)
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
      Pass_Both = CM_Pass & CZ_Pass,
      Special_Gene = orthogroup %in% special_gene_ids
    ) %>%
    dplyr::filter(Pass_Either | Special_Gene)
  
  cat("Number of retained effector orthogroups:", nrow(filtered_data), "\n")
  cat("CM-passed entries:", sum(filtered_data$CM_Pass, na.rm = TRUE), "\n")
  cat("CZ-passed entries:", sum(filtered_data$CZ_Pass, na.rm = TRUE), "\n")
  cat("Entries passing both CM and CZ criteria:",
      sum(filtered_data$Pass_Both, na.rm = TRUE), "\n")
  cat("Force-retained special genes:",
      sum(filtered_data$Special_Gene, na.rm = TRUE), "\n")
  
  return(filtered_data)
}


# ==============================================================================
# 10. Draw fire-style effector heatmap
# ==============================================================================

create_effector_heatmap_fire <- function(data,
                                         max_genes = 150,
                                         special_genes = NULL,
                                         cluster_rows = TRUE,
                                         output_dir = ".",
                                         fig_width = 10,
                                         fig_height = 8,
                                         tiff_res = 600) {
  cat("\nDrawing fire-style effector heatmap...\n")
  
  if (is.null(special_genes)) {
    special_genes <- character(0)
  }
  
  available_expr_cols <- grep(
    "^(CM|CZ)_(in vitro|0dpi|1dpi|3dpi|6dpi)$",
    colnames(data),
    value = TRUE
  )
  available_expr_cols <- sort_expression_columns(available_expr_cols)
  
  cat("Expression columns used in heatmap:",
      paste(available_expr_cols, collapse = ", "), "\n")
  
  if (length(available_expr_cols) == 0) {
    stop("No valid expression columns were found for heatmap.")
  }
  
  expr_data <- data %>%
    dplyr::filter(!is.na(orthogroup), orthogroup != "") %>%
    dplyr::rowwise() %>%
    dplyr::filter(
      sum(dplyr::c_across(dplyr::all_of(available_expr_cols)), na.rm = TRUE) > 0 |
        orthogroup %in% names(special_genes)
    ) %>%
    dplyr::ungroup()
  
  if (nrow(expr_data) == 0) {
    stop("No valid expression data are available for heatmap.")
  }
  
  special_data <- expr_data %>%
    dplyr::filter(orthogroup %in% names(special_genes))
  
  other_data <- expr_data %>%
    dplyr::filter(!orthogroup %in% names(special_genes))
  
  remaining_slots <- max_genes - nrow(special_data)
  
  if (remaining_slots > 0 && nrow(other_data) > 0) {
    other_data <- other_data %>%
      dplyr::mutate(
        Heatmap_Score = (
          (CM_3dpi + CM_6dpi + CZ_3dpi + CZ_6dpi) -
            (`CM_in vitro` + CM_0dpi + `CZ_in vitro` + CZ_0dpi)
        )
      ) %>%
      dplyr::arrange(dplyr::desc(Heatmap_Score)) %>%
      dplyr::slice_head(n = remaining_slots)
    
  } else if (remaining_slots <= 0) {
    other_data <- other_data[0, ]
  }
  
  final_expr_data <- dplyr::bind_rows(special_data, other_data) %>%
    dplyr::distinct(orthogroup, .keep_all = TRUE) %>%
    dplyr::mutate(
      Gene_Name = ifelse(
        orthogroup %in% names(special_genes),
        special_genes[orthogroup],
        ""
      ),
      Display_Label = make_display_label(orthogroup, special_genes)
    )
  
  cat(
    "Number of genes shown in heatmap:", nrow(final_expr_data),
    "(special genes:", nrow(special_data),
    "; other genes:", nrow(other_data), ")\n"
  )
  
  expr_matrix_data <- final_expr_data %>%
    dplyr::select(Display_Label, dplyr::all_of(available_expr_cols))
  
  expr_matrix <- convert_to_matrix_with_rownames(expr_matrix_data, "Display_Label")
  expr_matrix[is.na(expr_matrix)] <- 0
  expr_matrix[is.infinite(expr_matrix)] <- 0
  
  expr_matrix_log <- log2(expr_matrix + 1)
  
  row_annotation_df <- final_expr_data %>%
    dplyr::select(Display_Label, Presence_type) %>%
    as.data.frame()
  
  rownames(row_annotation_df) <- row_annotation_df$Display_Label
  row_annotation_df <- row_annotation_df[rownames(expr_matrix_log), , drop = FALSE]
  
  annotation_colors <- list(
    Presence_type = create_presence_type_colors(row_annotation_df$Presence_type)
  )
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Presence_type = row_annotation_df$Presence_type,
    col = annotation_colors,
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    annotation_name_rot = 45,
    simple_anno_size = grid::unit(0.45, "cm"),
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
    "0dpi" = "#CFE2F3",
    "1dpi" = "#FFF2CC",
    "3dpi" = "#FFE599",
    "6dpi" = "#FFD966"
  )
  
  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    Strain = col_info$Strain,
    Time_Point = col_info$Time_Point,
    col = list(
      Strain = c("CM" = "#FF6B6B", "CZ" = "#4ECDC4"),
      Time_Point = time_point_colors
    ),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
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
    row_names_gp = grid::gpar(fontsize = 7),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 45,
    left_annotation = row_ha,
    top_annotation = col_ha,
    column_title = NULL,
    width = grid::unit(12, "cm"),
    height = grid::unit(max(10, nrow(expr_matrix_log) * 0.22), "cm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.5)
  )
  
  pdf_file <- file.path(
    output_dir,
    "Effector_heatmap_fire_low_early_high_late.pdf"
  )
  
  tiff_file <- file.path(
    output_dir,
    "Effector_heatmap_fire_low_early_high_late.tiff"
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
  
  return(ht)
}


# ==============================================================================
# 11. Run analysis
# ==============================================================================

message_header("Data preprocessing")

processed_effector_data <- process_effector_expression_data(effector_data)

cat("Processed data dimensions:",
    paste(dim(processed_effector_data), collapse = " x "), "\n")

cat("\nPresence_type distribution:\n")
print(table(processed_effector_data$Presence_type, useNA = "ifany"))

cat("\nPreview of mean-expression columns:\n")
expr_preview_cols <- c("orthogroup", desired_expr_order)
expr_preview_cols <- expr_preview_cols[
  expr_preview_cols %in% colnames(processed_effector_data)
]
print(utils::head(processed_effector_data[, expr_preview_cols], 5))

special_genes_in_data <- names(special_genes)[
  names(special_genes) %in% processed_effector_data$orthogroup
]

cat("\nSpecial gene check:\n")
cat("Special genes present in data:",
    paste(special_genes_in_data, collapse = ", "), "\n")
cat("Special genes missing from data:",
    paste(setdiff(names(special_genes), special_genes_in_data), collapse = ", "), "\n")


message_header("Filtering target effector orthogroups")

filtered_effector_data <- filter_low_early_high_late_effectors(
  processed_effector_data,
  early_low_threshold = early_low_threshold,
  late_high_threshold = late_high_threshold,
  fold_change_threshold = fold_change_threshold,
  special_gene_ids = names(special_genes)
)

if (nrow(filtered_effector_data) == 0) {
  stop("No effector orthogroups passed the filtering criteria.")
}


message_header("Exporting filtered results")

filtered_output <- filtered_effector_data %>%
  dplyr::mutate(
    Gene_Name = ifelse(
      orthogroup %in% names(special_genes),
      special_genes[orthogroup],
      ""
    ),
    Display_Label = make_display_label(orthogroup, special_genes),
    Special_Gene = orthogroup %in% names(special_genes)
  ) %>%
  dplyr::arrange(
    dplyr::desc(Special_Gene),
    dplyr::desc(Pass_Both),
    dplyr::desc(CM_Late_vs_Early_FC + CZ_Late_vs_Early_FC)
  )

csv_file <- file.path(output_dir, "effector_filtered_low_early_high_late.csv")
xlsx_file <- file.path(output_dir, "effector_filtered_low_early_high_late.xlsx")

readr::write_csv(filtered_output, csv_file)
openxlsx::write.xlsx(filtered_output, xlsx_file, rowNames = FALSE)

cat("Filtered result files saved:\n")
cat(" -", csv_file, "\n")
cat(" -", xlsx_file, "\n")


message_header("Drawing heatmap")

heatmap_fire <- create_effector_heatmap_fire(
  filtered_effector_data,
  max_genes = max_genes_for_heatmap,
  special_genes = special_genes,
  cluster_rows = TRUE,
  output_dir = output_dir,
  fig_width = fig_width,
  fig_height = fig_height,
  tiff_res = tiff_res
)


# ==============================================================================
# 12. Save analysis report
# ==============================================================================

message_header("Saving analysis report")

report_file <- file.path(output_dir, "effector_heatmap_analysis_report.txt")

sink(report_file)

cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("Effector expression heatmap analysis report\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

cat("Analysis time:", as.character(Sys.time()), "\n")
cat("Input file:", input_file, "\n\n")

cat("1. Filtering criteria:\n")
cat(sprintf("   - Early low-expression threshold: mean(in vitro, 0 dpi) <= %s\n",
            early_low_threshold))
cat(sprintf("   - Late high-expression threshold: mean(3 dpi, 6 dpi) >= %s\n",
            late_high_threshold))
cat(sprintf("   - Late/Early fold-change threshold: >= %s\n",
            fold_change_threshold))

cat("\n2. Data summary:\n")
cat(sprintf("   - Total processed effector orthogroups: %d\n",
            nrow(processed_effector_data)))
cat(sprintf("   - Retained effector orthogroups: %d\n",
            nrow(filtered_effector_data)))
cat(sprintf("   - CM-passed entries: %d\n",
            sum(filtered_effector_data$CM_Pass, na.rm = TRUE)))
cat(sprintf("   - CZ-passed entries: %d\n",
            sum(filtered_effector_data$CZ_Pass, na.rm = TRUE)))
cat(sprintf("   - Entries passing both CM and CZ criteria: %d\n",
            sum(filtered_effector_data$Pass_Both, na.rm = TRUE)))
cat(sprintf("   - Force-retained special genes: %d\n",
            sum(filtered_effector_data$orthogroup %in% names(special_genes))))

cat("\n3. Presence_type distribution:\n")
presence_stats <- table(filtered_effector_data$Presence_type)

for (i in names(presence_stats)) {
  cat(sprintf(
    "   %s: %d entries (%.1f%%)\n",
    i,
    presence_stats[i],
    presence_stats[i] / nrow(filtered_effector_data) * 100
  ))
}

cat("\n4. Special gene retention:\n")
for (i in seq_along(special_genes)) {
  og <- names(special_genes)[i]
  gene_name <- special_genes[i]
  status <- ifelse(og %in% filtered_effector_data$orthogroup, "yes", "no")
  cat(sprintf("   %s | %s | retained: %s\n", og, gene_name, status))
}

cat("\n5. Heatmap column order:\n")
for (i in desired_expr_order) {
  cat(sprintf("   - %s\n", i))
}

cat("\n6. Output files:\n")
cat("   - Effector_heatmap_fire_low_early_high_late.pdf\n")
cat("   - Effector_heatmap_fire_low_early_high_late.tiff\n")
cat("   - effector_filtered_low_early_high_late.csv\n")
cat("   - effector_filtered_low_early_high_late.xlsx\n")
cat("   - effector_heatmap_analysis_report.txt\n")
cat("   - sessionInfo_effector_heatmap.txt\n")

cat("\n7. Heatmap notes:\n")
cat("   - Row annotation: Presence_type\n")
cat("   - Presence_type annotation title: 45 degrees, bold, not italicized\n")
cat("   - Ordinary genes are labeled as orthogroup IDs\n")
cat("   - Special genes are labeled as orthogroup | gene name\n")
cat("   - Top column annotations: Strain and Time_Point, not italicized\n")
cat("   - Color scheme: fire-style scale\n")
cat("   - Heatmap legend title: Log2(FPKM+1)\n")
cat("   - Heatmap title: removed\n")

cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

sink()

cat("Analysis report saved:\n")
cat(" -", report_file, "\n")


# ==============================================================================
# 13. Save session information
# ==============================================================================

session_file <- file.path(output_dir, "sessionInfo_effector_heatmap.txt")

sink(session_file)
cat("Session information\n")
cat("===================\n\n")
print(sessionInfo())
sink()

cat("Session information saved:\n")
cat(" -", session_file, "\n")


# ==============================================================================
# 14. Console summary
# ==============================================================================

message_header("Effector expression heatmap analysis completed")

cat("Filtering criteria:\n")
cat(sprintf(" - Early low-expression threshold: mean(in vitro, 0 dpi) <= %s\n",
            early_low_threshold))
cat(sprintf(" - Late high-expression threshold: mean(3 dpi, 6 dpi) >= %s\n",
            late_high_threshold))
cat(sprintf(" - Late/Early fold-change threshold: >= %s\n",
            fold_change_threshold))

cat("\nData summary:\n")
cat(sprintf(" - Total processed effector orthogroups: %d\n",
            nrow(processed_effector_data)))
cat(sprintf(" - Retained effector orthogroups: %d\n",
            nrow(filtered_effector_data)))
cat(sprintf(" - CM-passed entries: %d\n",
            sum(filtered_effector_data$CM_Pass, na.rm = TRUE)))
cat(sprintf(" - CZ-passed entries: %d\n",
            sum(filtered_effector_data$CZ_Pass, na.rm = TRUE)))
cat(sprintf(" - Entries passing both CM and CZ criteria: %d\n",
            sum(filtered_effector_data$Pass_Both, na.rm = TRUE)))
cat(sprintf(" - Force-retained special genes: %d\n",
            sum(filtered_effector_data$orthogroup %in% names(special_genes))))

cat("\nOutput directory:\n")
cat(" -", output_dir, "\n")

cat("\nOutput files:\n")
cat(" -", csv_file, "\n")
cat(" -", xlsx_file, "\n")
cat(" -", file.path(output_dir, "Effector_heatmap_fire_low_early_high_late.pdf"), "\n")
cat(" -", file.path(output_dir, "Effector_heatmap_fire_low_early_high_late.tiff"), "\n")
cat(" -", report_file, "\n")
cat(" -", session_file, "\n")

cat("\nAnalysis completed successfully.\n")
