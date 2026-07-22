#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_go_enrichment_infection_up_genes.R
#
# Purpose:
#   Perform GO enrichment analysis for genes upregulated during infection
#   relative to 0 dpi.
#
# Input files:
#   - data/background_GO.xlsx
#   - data/CM_infection_UP.csv
#   - data/CZ_infection_UP.csv
#
# Output:
#   - GO enrichment result tables
#   - Plot data tables
#   - Bubble plots
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
output_dir <- "results/go_infection_up"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

background_go_file <- file.path(input_dir, "background_GO.xlsx")

use_existing_plotdata <- FALSE

species_config <- list(
  CM = list(
    up_file      = file.path(input_dir, "CM_infection_UP.csv"),
    protein_col  = "CM_Gene.x",
    eggnog_col   = "CM_eggNOG_GO",
    interpro_col = "CM_InterPro_GO"
  ),
  CZ = list(
    up_file      = file.path(input_dir, "CZ_infection_UP.csv"),
    protein_col  = "CZ_Gene.x",
    eggnog_col   = "CZ_eggNOG_GO",
    interpro_col = "CZ_InterPro_GO"
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
    "proton-transporting ATP synthase complex" = "ATP synthase complex",
    "exosome (RNase complex)" = "exosome complex",
    "inner mitochondrial membrane protein complex" = "inner mito. membrane protein complex",
    "proton-transporting ATP synthase activity" = "ATP synthase activity",
    "ribonucleoprotein complex binding" = "RNP complex binding",
    "RNA-templated DNA biosynthetic process" = "RNA-templated DNA biosynthesis",
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

build_gene_go_map <- function(bg_df, protein_col, eggnog_col, interpro_col) {
  gene_to_go <- list()
  
  bg_sp <- bg_df[!is.na(bg_df[[protein_col]]) & bg_df[[protein_col]] != "", ]
  
  for (i in seq_len(nrow(bg_sp))) {
    gene_id <- bg_sp[[protein_col]][i]
    go_terms <- extract_go_terms(bg_sp[[eggnog_col]][i], bg_sp[[interpro_col]][i])
    
    if (length(go_terms) > 0) {
      gene_to_go[[gene_id]] <- unique(c(gene_to_go[[gene_id]], go_terms))
    }
  }
  
  gene_to_go
}

perform_go_enrichment <- function(target_genes, gene_go_map) {
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
    
    data.frame(
      GO_ID = go_id,
      Count = target_in_go,
      Target_total = n_target,
      Background_in_GO = length(genes_in_go),
      Background_total = n_bg,
      GeneRatio = target_in_go / n_target,
      P_value = fisher_result$p.value,
      Genes = paste(intersect(target_with_go, genes_in_go), collapse = "/"),
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

prepare_top_plot_data <- function(enrich_df, dataset_name, n_bp = 10, n_cc = 10, n_mf = 10) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  
  enrich_df <- enrich_df %>%
    arrange(Ontology, GO_term, P_value) %>%
    distinct(Ontology, GO_term, .keep_all = TRUE)
  
  top_df <- bind_rows(
    enrich_df %>% filter(Ontology == "BP") %>% arrange(P_value) %>% slice_head(n = n_bp),
    enrich_df %>% filter(Ontology == "CC") %>% arrange(P_value) %>% slice_head(n = n_cc),
    enrich_df %>% filter(Ontology == "MF") %>% arrange(P_value) %>% slice_head(n = n_mf)
  )
  
  if (nrow(top_df) == 0) return(NULL)
  
  top_df <- top_df %>%
    mutate(
      GO_term_plot = shorten_go_term(GO_term),
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
      plot.margin = margin(10, 20, 80, 70)
    )
}

# -----------------------------
# Run analysis
# -----------------------------

check_file(background_go_file)
bg_go <- readxl::read_excel(background_go_file)

for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  check_file(cfg$up_file)
  
  message("\nRunning infection-up GO enrichment for ", sp)
  
  up_df <- read.csv(cfg$up_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  if (!"Geneid" %in% colnames(up_df)) {
    stop("Column 'Geneid' was not found in ", cfg$up_file)
  }
  
  target_genes <- unique(up_df$Geneid)
  
  gene_to_go <- build_gene_go_map(
    bg_df = bg_go,
    protein_col = cfg$protein_col,
    eggnog_col = cfg$eggnog_col,
    interpro_col = cfg$interpro_col
  )
  
  enrich_raw <- perform_go_enrichment(target_genes, gene_to_go)
  enrich_final <- annotate_and_filter_go(enrich_raw)
  
  if (is.null(enrich_final) || nrow(enrich_final) == 0) {
    warning("No GO enrichment result for ", sp)
    next
  }
  
  out_all <- file.path(output_dir, paste0(sp, "_infection_UP_GO_all.csv"))
  write.csv(enrich_final, out_all, row.names = FALSE, quote = FALSE)
  
  plotdata_file <- file.path(output_dir, paste0(sp, "_infection_UP_GO_plotdata.csv"))
  
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
  
  ggsave(file.path(output_dir, paste0(sp, "_infection_UP_GO_top_plot.pdf")),
         plot = p, width = 16, height = 8.5, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(sp, "_infection_UP_GO_top_plot.png")),
         plot = p, width = 14, height = 10, dpi = 600, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(sp, "_infection_UP_GO_top_plot.tiff")),
         plot = p, width = 14, height = 10, dpi = 600, limitsize = FALSE)
}

sink(file.path(output_dir, "sessionInfo_go_infection_up.txt"))
print(sessionInfo())
sink()

message("\nInfection-up GO enrichment analysis completed.")
