#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_go_enrichment_specific_expressed_genes.R
#
# Purpose:
#   Perform GO enrichment analysis for CM-specific and CZ-specific expressed genes.
#
# Target gene definition:
#   species-specific genes intersected with genes expressed in at least one sample
#   using FPKM > expression_threshold.
#
# Input:
#   - data/background_GO.xlsx
#   - data/CM_specific.csv
#   - data/CZ_specific.csv
#   - data/CM_FPKM.csv
#   - data/CZ_FPKM.csv
#
# Output:
#   - GO enrichment result tables
#   - target gene to orthogroup mapping tables
#   - bubble plots
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(GO.db)
  library(dplyr)
  library(tidyr)
  library(tidytext)
})

# -----------------------------
# Parameters
# -----------------------------

input_dir <- "data"
output_dir <- "results/go_specific_expressed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

background_go_file <- file.path(input_dir, "background_GO.xlsx")
expression_threshold <- 0
use_existing_plotdata <- FALSE

species_config <- list(
  CM = list(
    specific_file = file.path(input_dir, "CM_specific.csv"),
    fpkm_file     = file.path(input_dir, "CM_FPKM.csv"),
    protein_col   = "CM_Gene.x",
    eggnog_col    = "CM_eggNOG_GO",
    interpro_col  = "CM_InterPro_GO"
  ),
  CZ = list(
    specific_file = file.path(input_dir, "CZ_specific.csv"),
    fpkm_file     = file.path(input_dir, "CZ_FPKM.csv"),
    protein_col   = "CZ_Gene.x",
    eggnog_col    = "CZ_eggNOG_GO",
    interpro_col  = "CZ_InterPro_GO"
  )
)

generic_go_terms <- c(
  "biological_process", "biological process", "metabolic process", "cellular process",
  "cellular metabolic process", "single-organism process", "single organism process",
  "biological regulation", "regulation of biological process",
  "cellular component organization", "cellular component organization or biogenesis",
  "localization", "establishment of localization", "transport",
  "response to stimulus", "signaling", "cell communication",
  "cellular_component", "cellular component", "cell", "cell part",
  "intracellular", "intracellular part", "organelle", "membrane", "membrane part",
  "organelle part", "molecular_function", "molecular function", "binding",
  "catalytic activity", "protein binding", "ion binding", "nucleotide binding",
  "nucleic acid binding", "hydrolase activity", "transferase activity",
  "transporter activity", "Unknown", "unknown", "obsolete", "..."
)

# -----------------------------
# Helper functions
# -----------------------------

check_file <- function(file) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }
}

extract_go_terms <- function(eggnog_go, interpro_go) {
  go_terms <- character(0)
  
  if (!is.na(eggnog_go) && eggnog_go != "") {
    go_terms <- c(go_terms, stringr::str_extract_all(eggnog_go, "GO:\\d+")[[1]])
  }
  
  if (!is.na(interpro_go) && interpro_go != "") {
    go_terms <- c(go_terms, stringr::str_extract_all(interpro_go, "GO:\\d+")[[1]])
  }
  
  unique(go_terms)
}

go_id_to_term <- function(go_ids) {
  sapply(go_ids, function(go_id) {
    tryCatch({
      term <- Term(GOTERM[[go_id]])
      if (is.null(term) || is.na(term)) go_id else term
    }, error = function(e) go_id)
  })
}

go_id_to_ontology <- function(go_ids) {
  sapply(go_ids, function(go_id) {
    tryCatch({
      ont <- Ontology(GOTERM[[go_id]])
      if (is.null(ont) || is.na(ont)) "Unknown" else ont
    }, error = function(e) "Unknown")
  })
}

shorten_go_term <- function(term) {
  term_map <- c(
    "urea metabolic process" = "urea metabolism",
    "hydrogen peroxide catabolic process" = "hydrogen peroxide catabolism",
    "regulation of translational elongation" = "translational elongation",
    "nitrogen cycle metabolic process" = "nitrogen cycle metabolism",
    "proton motive force-driven ATP synthesis" = "ATP synthesis (PMF-driven)",
    "regulation of translational initiation" = "translational initiation",
    "Golgi to vacuole transport" = "Golgi-to-vacuole transport",
    "glutathione metabolic process" = "glutathione metabolism",
    "regulation of cytoplasmic translation" = "cytoplasmic translation",
    "regulation of cytoplasmic translational elongation" = "reg. cytoplasmic translational elongation",
    "proton-transporting ATP synthase complex" = "ATP synthase complex",
    "exosome (RNase complex)" = "exosome complex",
    "inner mitochondrial membrane protein complex" = "inner mito. membrane protein complex",
    "proton-transporting ATP synthase activity" = "ATP synthase activity",
    "ribonucleoprotein complex binding" = "RNP complex binding",
    "RNA-templated DNA biosynthetic process" = "RNA-templated DNA biosynthesis",
    "double-strand break repair via synthesis-dependent strand annealing" = "DSB repair via SDSA",
    "double-strand break repair via nonhomologous end joining" = "DSB repair via NHEJ",
    "double-strand break repair" = "DSB repair",
    "site of DNA damage" = "DNA damage site",
    "site of double-strand break" = "DSB site",
    "single-stranded DNA binding" = "ssDNA binding",
    "double-stranded DNA binding" = "dsDNA binding",
    "nucleobase-containing compound kinase activity" = "nucleobase compound kinase activity",
    "ATP-dependent activity, acting on DNA" = "ATP-dependent activity on DNA",
    "catalytic activity, acting on DNA" = "catalytic activity on DNA"
  )
  
  ifelse(term %in% names(term_map), term_map[term], term)
}

truncate_term <- function(term, max_length = 40) {
  ifelse(nchar(term) > max_length, paste0(substr(term, 1, max_length), "..."), term)
}

get_expressed_genes <- function(fpkm_file, gene_col = 1, threshold = 0) {
  check_file(fpkm_file)
  
  fpkm_df <- read.csv(fpkm_file, check.names = FALSE, stringsAsFactors = FALSE)
  gene_ids <- fpkm_df[[gene_col]]
  
  expr_cols <- setdiff(colnames(fpkm_df), colnames(fpkm_df)[gene_col])
  expr_mat <- fpkm_df[, expr_cols, drop = FALSE]
  expr_mat[] <- lapply(expr_mat, function(x) suppressWarnings(as.numeric(x)))
  
  max_expr <- apply(expr_mat, 1, function(x) max(x, na.rm = TRUE))
  expressed_genes <- gene_ids[max_expr > threshold]
  
  unique(expressed_genes[!is.na(expressed_genes) & expressed_genes != ""])
}

build_gene_maps <- function(bg_df, protein_col, eggnog_col, interpro_col, orthogroup_col = "orthogroup") {
  gene_to_go <- list()
  gene_to_og <- list()
  
  bg_sp <- bg_df[!is.na(bg_df[[protein_col]]) & bg_df[[protein_col]] != "", ]
  
  for (i in seq_len(nrow(bg_sp))) {
    gene_id <- bg_sp[[protein_col]][i]
    og_id <- if (orthogroup_col %in% colnames(bg_sp)) bg_sp[[orthogroup_col]][i] else NA
    
    go_terms <- extract_go_terms(bg_sp[[eggnog_col]][i], bg_sp[[interpro_col]][i])
    
    if (length(go_terms) > 0) {
      gene_to_go[[gene_id]] <- unique(c(gene_to_go[[gene_id]], go_terms))
    }
    
    gene_to_og[[gene_id]] <- unique(c(gene_to_og[[gene_id]], og_id))
  }
  
  list(gene_to_go = gene_to_go, gene_to_og = gene_to_og)
}

perform_go_enrichment <- function(target_genes, gene_go_map, gene_og_map) {
  go_to_genes <- list()
  
  for (gene in names(gene_go_map)) {
    for (go_id in gene_go_map[[gene]]) {
      go_to_genes[[go_id]] <- unique(c(go_to_genes[[go_id]], gene))
    }
  }
  
  target_with_go <- intersect(target_genes, names(gene_go_map))
  background_with_go <- names(gene_go_map)
  
  n_target <- length(target_with_go)
  n_bg <- length(background_with_go)
  
  message("  Target genes with GO: ", n_target, "/", length(target_genes))
  message("  Background genes with GO: ", n_bg)
  
  enrichment_list <- lapply(names(go_to_genes), function(go_id) {
    genes_in_go <- go_to_genes[[go_id]]
    
    target_in_go <- length(intersect(target_with_go, genes_in_go))
    target_not_in_go <- n_target - target_in_go
    bg_in_go_not_target <- length(setdiff(genes_in_go, target_with_go))
    bg_not_in_go_not_target <- n_bg - length(genes_in_go) - target_not_in_go
    
    if (target_in_go < 2 || length(genes_in_go) < 5) return(NULL)
    
    contingency_table <- matrix(
      c(target_in_go, target_not_in_go, bg_in_go_not_target, bg_not_in_go_not_target),
      nrow = 2,
      byrow = TRUE
    )
    
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    hit_genes <- intersect(target_with_go, genes_in_go)
    
    hit_ogs <- unique(unlist(gene_og_map[hit_genes]))
    hit_ogs <- hit_ogs[!is.na(hit_ogs) & hit_ogs != ""]
    
    data.frame(
      GO_ID = go_id,
      Count = target_in_go,
      Target_total = n_target,
      Background_in_GO = length(genes_in_go),
      Background_total = n_bg,
      GeneRatio = target_in_go / n_target,
      P_value = fisher_result$p.value,
      Genes = paste(hit_genes, collapse = "/"),
      Orthogroups = paste(hit_ogs, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  enrichment_list <- enrichment_list[!sapply(enrichment_list, is.null)]
  if (length(enrichment_list) == 0) return(NULL)
  
  enrichment_df <- do.call(rbind, enrichment_list)
  enrichment_df$p.adjusted <- p.adjust(enrichment_df$P_value, method = "BH")
  enrichment_df <- enrichment_df[order(enrichment_df$P_value), ]
  
  enrichment_df
}

annotate_and_filter_go <- function(enrich_df) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  enrich_df$GO_term <- as.character(go_id_to_term(enrich_df$GO_ID))
  enrich_df$Ontology <- as.character(go_id_to_ontology(enrich_df$GO_ID))
  
  enrich_df <- enrich_df[
    !tolower(enrich_df$GO_term) %in% tolower(generic_go_terms) &
      enrich_df$Ontology != "Unknown",
  ]
  
  if (nrow(enrich_df) == 0) return(NULL)
  enrich_df
}

prepare_top_plot_data <- function(enrich_df, species_name) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  enrich_df <- enrich_df %>%
    arrange(Ontology, GO_term, P_value) %>%
    distinct(Ontology, GO_term, .keep_all = TRUE)
  
  if (species_name == "CM") {
    n_bp <- 8
    n_cc <- 6
    n_mf <- 6
  } else {
    n_bp <- 10
    n_cc <- 10
    n_mf <- 10
  }
  
  top_df <- bind_rows(
    enrich_df %>% filter(Ontology == "BP") %>% arrange(P_value) %>% slice_head(n = n_bp),
    enrich_df %>% filter(Ontology == "CC") %>% arrange(P_value) %>% slice_head(n = n_cc),
    enrich_df %>% filter(Ontology == "MF") %>% arrange(P_value) %>% slice_head(n = n_mf)
  )
  
  if (nrow(top_df) == 0) return(NULL)
  
  top_df <- top_df %>%
    mutate(
      GO_term_plot = truncate_term(shorten_go_term(GO_term), max_length = 40),
      Ontology_label = factor(Ontology, levels = c("BP", "CC", "MF")),
      logP = -log10(P_value)
    )
  
  duplicated_terms <- duplicated(top_df$GO_term_plot) | duplicated(top_df$GO_term_plot, fromLast = TRUE)
  top_df$GO_term_plot[duplicated_terms] <- paste0(
    top_df$GO_term_plot[duplicated_terms],
    " (",
    top_df$GO_ID[duplicated_terms],
    ")"
  )
  
  top_df$GO_term_plot <- tidytext::reorder_within(
    top_df$GO_term_plot,
    -top_df$logP,
    top_df$Ontology_label
  )
  
  top_df
}

create_bubble_plot <- function(plot_df) {
  size_breaks <- sort(unique(plot_df$Count))
  if (length(size_breaks) > 5) {
    size_breaks <- pretty(range(plot_df$Count), n = 4)
  }
  
  y_max <- max(plot_df$GeneRatio, na.rm = TRUE)
  y_upper <- ifelse(y_max < 0.15, y_max * 1.35, y_max * 1.12)
  
  ggplot(plot_df, aes(x = GO_term_plot, y = GeneRatio)) +
    geom_point(aes(size = Count, color = -log10(P_value)), alpha = 0.95) +
    facet_grid(. ~ Ontology_label, scales = "free_x", space = "free_x") +
    tidytext::scale_x_reordered(expand = expansion(add = 0.6)) +
    scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0.01, 0.02))) +
    scale_colour_gradient(low = "#63B8FF", high = "#F94144", name = expression(-log[10](pvalue))) +
    scale_size_continuous(name = "Gene counts", range = c(5, 16), breaks = size_breaks) +
    labs(x = NULL, y = "Gene ratio") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
      strip.background = element_rect(fill = "grey90", color = "grey60"),
      strip.text.x = element_text(size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.05, size = 14),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      plot.title = element_blank(),
      plot.margin = margin(10, 20, 60, 130)
    )
}

# -----------------------------
# Run analysis
# -----------------------------

check_file(background_go_file)
bg_go <- readxl::read_excel(background_go_file)

for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  
  check_file(cfg$specific_file)
  check_file(cfg$fpkm_file)
  
  message("\nRunning specific-expressed GO enrichment for ", sp)
  
  specific_df <- read.csv(cfg$specific_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  if (!"gene_id" %in% colnames(specific_df)) {
    stop("Column 'gene_id' was not found in ", cfg$specific_file)
  }
  
  specific_genes <- unique(specific_df$gene_id)
  expressed_genes <- get_expressed_genes(cfg$fpkm_file, gene_col = 1, threshold = expression_threshold)
  target_genes <- intersect(specific_genes, expressed_genes)
  
  message("  Specific genes: ", length(specific_genes))
  message("  Expressed genes: ", length(expressed_genes))
  message("  Target genes: ", length(target_genes))
  
  maps <- build_gene_maps(
    bg_df = bg_go,
    protein_col = cfg$protein_col,
    eggnog_col = cfg$eggnog_col,
    interpro_col = cfg$interpro_col,
    orthogroup_col = "orthogroup"
  )
  
  enrich_raw <- perform_go_enrichment(target_genes, maps$gene_to_go, maps$gene_to_og)
  enrich_final <- annotate_and_filter_go(enrich_raw)
  
  if (is.null(enrich_final) || nrow(enrich_final) == 0) {
    warning("No GO enrichment result for ", sp)
    next
  }
  
  out_all <- file.path(output_dir, paste0(sp, "_specific_expressed_GO_all.csv"))
  write.csv(enrich_final, out_all, row.names = FALSE, quote = FALSE)
  
  target_og_df <- data.frame(
    GeneID = target_genes,
    Orthogroup = sapply(target_genes, function(g) {
      og <- unique(unlist(maps$gene_to_og[g]))
      og <- og[!is.na(og) & og != ""]
      paste(og, collapse = "/")
    }),
    stringsAsFactors = FALSE
  )
  
  out_target <- file.path(output_dir, paste0(sp, "_specific_expressed_genes_with_orthogroup.csv"))
  write.csv(target_og_df, out_target, row.names = FALSE, quote = FALSE)
  
  plotdata_file <- file.path(output_dir, paste0(sp, "_specific_expressed_GO_plotdata.csv"))
  
  if (use_existing_plotdata && file.exists(plotdata_file)) {
    plot_df <- read.csv(plotdata_file, check.names = FALSE, stringsAsFactors = FALSE)
    plot_df$Ontology_label <- factor(plot_df$Ontology_label, levels = c("BP", "CC", "MF"))
    plot_df$logP <- -log10(plot_df$P_value)
    plot_df$GO_term_plot <- tidytext::reorder_within(
      plot_df$GO_term_plot,
      -plot_df$logP,
      plot_df$Ontology_label
    )
  } else {
    plot_df <- prepare_top_plot_data(enrich_final, sp)
    write.csv(plot_df, plotdata_file, row.names = FALSE, quote = FALSE)
  }
  
  p <- create_bubble_plot(plot_df)
  
  ggsave(file.path(output_dir, paste0(sp, "_specific_expressed_GO_top_plot.pdf")),
         plot = p, width = 16, height = 8, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(sp, "_specific_expressed_GO_top_plot.png")),
         plot = p, width = 16, height = 8, dpi = 600, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(sp, "_specific_expressed_GO_top_plot.tiff")),
         plot = p, width = 16, height = 8, dpi = 600, limitsize = FALSE)
}

sink(file.path(output_dir, "sessionInfo_go_specific_expressed.txt"))
print(sessionInfo())
sink()

message("\nSpecific-expressed GO enrichment analysis completed.")
