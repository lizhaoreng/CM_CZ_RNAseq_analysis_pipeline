#!/usr/bin/env Rscript

# ==============================================================================
# 05_analyze_dnds_functional_categories.R
#
# Analyze evolutionary rates of effector candidates, CWDEs, and background genes.
#
# Inputs:
#   - results/dnds_enhanced_clean.csv
#   - input/effector_one_to_one_list.txt
#   - input/CWDE_one_to_one_list.txt
#
# Outputs:
#   - Summary statistics
#   - Wilcoxon test results
#   - Fisher enrichment results
#   - Multivariable regression results
#   - Publication-ready figures
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggpubr)
  library(patchwork)
  library(rstatix)
  library(broom)
  library(scales)
  library(Cairo)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.yaml")

config <- yaml::read_yaml(config_file)

work_dir <- config$output$work_dir
figure_dir <- config$output$figures_dir

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

file_dnds <- file.path(work_dir, "dnds_enhanced_clean.csv")
file_eff <- config$input$effector_list
file_cwde <- config$input$cwde_list

omega_max_cutoff <- config$quality_filter$max_omega
ds_max_cutoff <- config$quality_filter$max_dS
alpha_level <- config$statistics$alpha
fdr_method <- config$statistics$fdr_method

cat("============================================================\n")
cat("Functional category dN/dS analysis\n")
cat("============================================================\n")

if (!file.exists(file_dnds)) stop("Missing dN/dS file: ", file_dnds)
if (!file.exists(file_eff)) stop("Missing effector list: ", file_eff)
if (!file.exists(file_cwde)) stop("Missing CWDE list: ", file_cwde)

dnds <- data.table::fread(file_dnds)
eff_list <- data.table::fread(file_eff, header = TRUE, sep = "\t")
cwde_list <- data.table::fread(file_cwde, header = TRUE, sep = "\t")

colnames(eff_list)[1] <- "Orthogroup"
colnames(cwde_list)[1] <- "Orthogroup"

eff_ogs <- unique(eff_list$Orthogroup)
cwde_ogs <- unique(cwde_list$Orthogroup)

cat("Input dN/dS records:", nrow(dnds), "\n")
cat("Effector orthogroups:", length(eff_ogs), "\n")
cat("CWDE orthogroups:", length(cwde_ogs), "\n")

d <- dnds %>%
  mutate(
    Category = case_when(
      Orthogroup %in% eff_ogs ~ "Effector",
      Orthogroup %in% cwde_ogs ~ "CWDE",
      TRUE ~ "Background"
    ),
    Category = factor(Category, levels = c("Background", "Effector", "CWDE")),
    omega_num = as.numeric(omega),
    dN_num = as.numeric(dN),
    dS_num = as.numeric(dS),
    Length = as.numeric(Effective_Codons),
    GapPct = as.numeric(Gap_Percentage),
    GC = as.numeric(Avg_GC_Content),
    IDY = as.numeric(Sequence_Identity)
  ) %>%
  filter(
    QC_Pass == TRUE,
    !is.na(omega_num),
    !is.na(dS_num),
    dS_num > 0,
    dS_num <= ds_max_cutoff,
    omega_num > 0,
    omega_num < omega_max_cutoff
  )

# Convert sequence identity from 0-1 scale to percent if needed
if (max(d$IDY, na.rm = TRUE) <= 1.2) {
  d <- d %>% mutate(IDY = IDY * 100)
}

cat("\nFinal category distribution:\n")
print(table(d$Category))

# -----------------------------
# Summary statistics
# -----------------------------

summary_stats <- d %>%
  group_by(Category) %>%
  summarise(
    n = n(),
    omega_median = median(omega_num, na.rm = TRUE),
    omega_mean = mean(omega_num, na.rm = TRUE),
    omega_sd = sd(omega_num, na.rm = TRUE),
    omega_q25 = quantile(omega_num, 0.25, na.rm = TRUE),
    omega_q75 = quantile(omega_num, 0.75, na.rm = TRUE),
    omega_min = min(omega_num, na.rm = TRUE),
    omega_max = max(omega_num, na.rm = TRUE),
    .groups = "drop"
  )

data.table::fwrite(
  summary_stats,
  file.path(figure_dir, "Table1_dnds_summary_statistics.tsv"),
  sep = "\t"
)

# -----------------------------
# Wilcoxon tests
# -----------------------------

comparison_pairs <- list(
  c("Background", "Effector"),
  c("Background", "CWDE"),
  c("Effector", "CWDE")
)

wilcox_results <- purrr::map_dfr(comparison_pairs, function(pair) {
  tmp <- d %>% filter(Category %in% pair)
  
  if (length(unique(tmp$Category)) != 2) {
    return(NULL)
  }
  
  test <- wilcox.test(omega_num ~ Category, data = tmp)
  
  group1 <- tmp$omega_num[tmp$Category == pair[1]]
  group2 <- tmp$omega_num[tmp$Category == pair[2]]
  
  n1 <- length(group1)
  n2 <- length(group2)
  
  dominance <- sum(outer(group2, group1, ">")) +
    0.5 * sum(outer(group2, group1, "=="))
  
  cliff_delta <- (2 * dominance) / (n1 * n2) - 1
  
  tibble(
    Comparison = paste(pair[2], "vs", pair[1]),
    Group1 = pair[1],
    Group2 = pair[2],
    n1 = n1,
    n2 = n2,
    median1 = median(group1, na.rm = TRUE),
    median2 = median(group2, na.rm = TRUE),
    W_statistic = as.numeric(test$statistic),
    p.value = test$p.value,
    cliff_delta = cliff_delta
  )
}) %>%
  mutate(
    q.value = p.adjust(p.value, method = fdr_method),
    significance = case_when(
      q.value < 0.001 ~ "***",
      q.value < 0.01 ~ "**",
      q.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

data.table::fwrite(
  wilcox_results,
  file.path(figure_dir, "Table2_wilcoxon_tests.tsv"),
  sep = "\t"
)

# -----------------------------
# Fisher enrichment tests
# -----------------------------

omega_top10_cut <- quantile(d$omega_num, 0.90, na.rm = TRUE)
omega_top5_cut <- quantile(d$omega_num, 0.95, na.rm = TRUE)

run_fisher_enrichment <- function(threshold_name, threshold_function) {
  purrr::map_dfr(list(c("Background", "Effector"), c("Background", "CWDE")), function(pair) {
    tmp <- d %>%
      filter(Category %in% pair) %>%
      mutate(high_omega = threshold_function(omega_num))
    
    contingency <- table(tmp$Category, tmp$high_omega)
    
    if (nrow(contingency) != 2 || ncol(contingency) != 2) {
      return(NULL)
    }
    
    fisher_result <- fisher.test(contingency)
    props <- prop.table(contingency, 1)
    
    tibble(
      Comparison = paste(pair[2], "vs", pair[1]),
      Threshold = threshold_name,
      n_background = sum(tmp$Category == pair[1]),
      n_functional = sum(tmp$Category == pair[2]),
      prop_background_high = props[pair[1], "TRUE"],
      prop_functional_high = props[pair[2], "TRUE"],
      odds_ratio = as.numeric(fisher_result$estimate),
      odds_ratio_low = fisher_result$conf.int[1],
      odds_ratio_high = fisher_result$conf.int[2],
      p.value = fisher_result$p.value
    )
  })
}

fisher_results <- bind_rows(
  run_fisher_enrichment("Top10%", function(x) x >= omega_top10_cut),
  run_fisher_enrichment("Top5%", function(x) x >= omega_top5_cut),
  run_fisher_enrichment("omega>0.5", function(x) x > 0.5),
  run_fisher_enrichment("omega>1", function(x) x > 1)
) %>%
  mutate(
    q.value = p.adjust(p.value, method = fdr_method),
    significance = case_when(
      q.value < 0.001 ~ "***",
      q.value < 0.01 ~ "**",
      q.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

data.table::fwrite(
  fisher_results,
  file.path(figure_dir, "Table3_fisher_enrichment_tests.tsv"),
  sep = "\t"
)

# -----------------------------
# Multivariable regression
# -----------------------------

d_reg <- d %>%
  filter(
    is.finite(omega_num),
    omega_num > 0,
    is.finite(Length),
    is.finite(IDY),
    is.finite(GC),
    is.finite(GapPct)
  ) %>%
  mutate(
    log_omega = log10(omega_num),
    log_length = log10(Length + 1)
  )

fit <- lm(log_omega ~ Category + log_length + IDY + GC + GapPct, data = d_reg)

regression_results <- broom::tidy(fit, conf.int = TRUE) %>%
  mutate(
    term_clean = case_when(
      term == "CategoryEffector" ~ "Effector vs Background",
      term == "CategoryCWDE" ~ "CWDE vs Background",
      TRUE ~ term
    )
  )

data.table::fwrite(
  regression_results,
  file.path(figure_dir, "Table4_multivariable_regression.tsv"),
  sep = "\t"
)

# -----------------------------
# Figures
# -----------------------------

theme_journal <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold")
  )

category_colors <- c(
  "Background" = "#2166AC",
  "Effector" = "#D6604D",
  "CWDE" = "#5AAE61"
)

category_shapes <- c(
  "Background" = 16,
  "Effector" = 17,
  "CWDE" = 15
)

p_density <- ggplot(d, aes(x = omega_num, color = Category, fill = Category)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  geom_vline(xintercept = omega_top10_cut, linetype = "longdash", color = "grey30") +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2),
    labels = c("0.01", "0.05", "0.1", "0.5", "1", "2")
  ) +
  scale_color_manual(values = category_colors) +
  scale_fill_manual(values = category_colors) +
  labs(
    x = expression(paste("dN/dS ratio (", omega, ", log scale)")),
    y = "Density",
    title = "Distribution of evolutionary rates"
  ) +
  theme_journal

p_box <- ggplot(d, aes(x = Category, y = omega_num, fill = Category, color = Category)) +
  geom_violin(alpha = 0.2, trim = TRUE, linewidth = 0.5) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8, linewidth = 0.5) +
  geom_jitter(aes(shape = Category), width = 0.12, alpha = 0.25, size = 0.8) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2),
    labels = c("0.01", "0.05", "0.1", "0.5", "1", "2")
  ) +
  scale_fill_manual(values = category_colors) +
  scale_color_manual(values = category_colors) +
  scale_shape_manual(values = category_shapes) +
  labs(
    x = NULL,
    y = expression(paste("dN/dS ratio (", omega, ", log scale)")),
    title = "Evolutionary rate comparison"
  ) +
  theme_journal +
  theme(legend.position = "none")

p_length <- ggplot(d, aes(x = Length, y = omega_num, color = Category, shape = Category)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, alpha = 0.2) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = category_colors) +
  scale_shape_manual(values = category_shapes) +
  labs(
    x = "Effective codons",
    y = expression(paste("dN/dS ratio (", omega, ")")),
    title = "Gene length vs evolutionary rate"
  ) +
  theme_journal

p_identity <- ggplot(d, aes(x = IDY, y = omega_num, color = Category, shape = Category)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, alpha = 0.2) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = category_colors) +
  scale_shape_manual(values = category_shapes) +
  labs(
    x = "Sequence identity (%)",
    y = expression(paste("dN/dS ratio (", omega, ")")),
    title = "Sequence identity vs evolutionary rate"
  ) +
  theme_journal

main_figure <- (p_density | p_box) / (p_length | p_identity)

ggsave(
  file.path(figure_dir, "dnds_functional_category_analysis.pdf"),
  main_figure,
  width = 14,
  height = 10,
  device = cairo_pdf
)

ggsave(
  file.path(figure_dir, "dnds_functional_category_analysis.png"),
  main_figure,
  width = 14,
  height = 10,
  dpi = 600
)

# Regression forest plot
reg_plot_data <- regression_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = factor(
      term_clean,
      levels = rev(c("Effector vs Background", "CWDE vs Background",
                     "log_length", "IDY", "GC", "GapPct"))
    ),
    significant = p.value < alpha_level
  )

p_regression <- ggplot(reg_plot_data, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, linewidth = 0.8) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey60")) +
  labs(
    x = expression(paste("Regression coefficient for log"[10], "(", omega, ")")),
    y = NULL,
    title = "Multivariable regression analysis"
  ) +
  theme_journal +
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir, "multivariable_regression_forest_plot.pdf"),
  p_regression,
  width = 9,
  height = 5,
  device = cairo_pdf
)

# Save processed data and session info
data.table::fwrite(
  d,
  file.path(figure_dir, "dnds_functional_category_data.tsv"),
  sep = "\t"
)

sink(file.path(figure_dir, "sessionInfo_dnds_functional_analysis.txt"))
print(sessionInfo())
sink()

cat("\nAnalysis completed successfully.\n")
cat("Output directory:", figure_dir, "\n")
