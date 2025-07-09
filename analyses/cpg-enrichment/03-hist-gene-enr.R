# Plot PBTA CPG PLP enrichment by plot group
# Ryan Corbett
# 2024

# Load libraries

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(scales)
library(ggsci)
library(ggpubr)
library(tidytext)

# Set directory path

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "cpg-enrichment")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))
source(file.path(analysis_dir, "util", "enrichment_functions.R"))

# Set file file paths

cpg_enr_gnomad_file <- file.path(input_dir, 
                                 "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_gene_gnomad_enrichment.tsv")
cpg_enr_pmbb_file <- file.path(input_dir, 
                               "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_gene_pmbb_enrichment.tsv")

opc_hist_file <- file.path(data_dir, 
                           "histologies.tsv")

cbtn_histologies_file <- file.path(root_dir, "analyses", 
                                   "collapse-tumor-histologies", "results", 
                                   "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

plp_all_exome_file <- file.path(root_dir, "analyses",
                                "bed-intersect", "results", 
                                "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded.tsv")


# Read in hist, plp, and enrichment by CPG files

cpgs <- read_lines(file.path(root_dir, "analyses", 
                             "oncokb-annotation",
                             "input", "cpg.txt"))

cpg_enr_gnomad <- read_tsv(cpg_enr_gnomad_file) 

cpg_enr_pmbb <- read_tsv(cpg_enr_pmbb_file) 

# define pineos as a separate group for plotting
hist <- read_tsv(cbtn_histologies_file) %>%
  dplyr::mutate(plot_group = case_when(
    grepl("Pineoblastoma", pathology_diagnosis) ~ "Pineoblastoma",
    TRUE ~ plot_group
  ))

## Calculate gnomad enrichment

#obtain sample counts per plot group 
hist_cts <- hist %>%
  count(plot_group) %>%
  dplyr::rename(total_cohort_size_case = n)

# Load autogvp output 
plp <- read_tsv(plp_all_exome_file) %>%
  dplyr::filter(gene_symbol_vep %in% cpgs) %>%
  distinct(Kids_First_Biospecimen_ID_normal, gene_symbol_vep, .keep_all = TRUE) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_normal, plot_group))

# Obtain P-LP carrier count by CPG and plot group
hist_plp_ct <- plp %>%
  count(gene_symbol_vep, plot_group) %>%
  dplyr::filter(n > 1) %>%
  dplyr::rename(count_with_plp_case = n) %>%
  right_join(hist_cts) %>%
  dplyr::mutate(count_without_plp_case = total_cohort_size_case - count_with_plp_case)
  
# Calculate plot-group level P-LP carrier enrichment relative to gnomAD for each CPG 
hist_cpg_enr_gnomad <- cpg_enr_gnomad %>%
  dplyr::filter(!is.na(gene_symbol_vep)) %>%
  dplyr::select(gene_symbol_vep, count_with_plp_control,
                total_cohort_size_control, count_without_plp_control) %>%
  left_join(hist_plp_ct) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  dplyr::mutate(p = NA_integer_,
                OR = NA_integer_,
                ci.int1 = NA_integer_,
                ci.int2 = NA_integer_)

hist_cpg_enr_gnomad <- calculate_enrichment(hist_cpg_enr_gnomad)

# multiple test correction 
hist_cpg_enr_gnomad <- hist_cpg_enr_gnomad %>%
  group_by(plot_group) %>%
  dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))

# Repeat enrichment calculations for PMBB controls

# calculate PBTA enrichment for P-LP variants by gene and plot group from PMBB controls
hist_cpg_enr_pmbb <- cpg_enr_pmbb %>%
  dplyr::filter(!is.na(gene_symbol_vep)) %>%
  dplyr::select(gene_symbol_vep, count_with_plp_control,
                total_cohort_size_control, count_without_plp_control) %>%
  left_join(hist_plp_ct) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  dplyr::mutate(p = NA_integer_,
                OR = NA_integer_,
                ci.int1 = NA_integer_,
                ci.int2 = NA_integer_)

hist_cpg_enr_pmbb <- calculate_enrichment(hist_cpg_enr_pmbb)

# mulitiple test correction
hist_cpg_enr_pmbb <- hist_cpg_enr_pmbb %>%
  group_by(plot_group) %>%
  dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))

# Create 'dummy' data frame with PBTA enrichment values to use as reference in plots
hist_cpg_enr_pbta <- hist_cpg_enr_pmbb %>%
  dplyr::select(gene_symbol_vep, plot_group, count_with_plp_case, count_without_plp_case, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PBTA",
         OR = NA_integer_,
         p = NA_integer_,
         ci.int1 = NA_integer_, 
         ci.int2 = NA_integer_,
         padj = NA_integer_) %>%
  dplyr::rename("n_plp" = "count_with_plp_case", 
                "n_no_plp" = "count_without_plp_case")

# format gnomAD and PMBB enrichment results for merging
hist_cpg_enr_gnomad <- hist_cpg_enr_gnomad %>%
  dplyr::select(gene_symbol_vep, plot_group, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "gnomAD") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

hist_cpg_enr_pmbb <- hist_cpg_enr_pmbb %>%
  dplyr::select(gene_symbol_vep, plot_group, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PMBB") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

# merge enrichment results and calculate percent of patients with CPG PLP
hist_cpg_enr_all <- hist_cpg_enr_pbta %>%
  bind_rows(hist_cpg_enr_gnomad, hist_cpg_enr_pmbb) %>%
  mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
         fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
  arrange(desc(-log10(p))) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
         hist_gene = glue::glue("{plot_group}: {gene_symbol_vep}"),
         hist_gene = factor(hist_gene, unique(hist_gene)),)

write_tsv(hist_cpg_enr_all,
          file.path(results_dir, "hist-gene-plp-enr-pbta-vs-pmbb-gnomad.tsv"))

# pull significnalty enriched CPGs relative to gnomAD and PMBB cohorts, and obtain those common to both sets
sig_hist_cpgs_gnomad <- hist_cpg_enr_all %>%
  filter(cohort == "gnomAD" & padj < 0.05 & padj > 0) %>%
  pull(hist_gene)

sig_hist_cpgs_pmbb <- hist_cpg_enr_all %>%
  filter(cohort == "PMBB" & padj < 0.05 & padj > 0) %>%
  pull(hist_gene)

sig_hist_cpgs_both <- intersect(sig_hist_cpgs_gnomad, sig_hist_cpgs_pmbb)

hist_cpg_enr_all <- hist_cpg_enr_all %>%
  dplyr::filter(hist_gene %in% sig_hist_cpgs_both)

# Loop through plot groups to plot sig results
plot_groups <- hist_cpg_enr_all %>%
  filter(hist_gene %in% sig_hist_cpgs_both) %>%
  pull(plot_group) %>%
  unique()

merged_plot <- list()

for (group in sort(plot_groups)) {
  
  # get number of significantly enriched genes in plot group to determine plot height
  n_sig <- sum(grepl(group, sig_hist_cpgs_both))
  
  # Create CPG enrichment FDR plot  
  
  pval_plot <- hist_cpg_enr_all %>%
    dplyr::filter(plot_group == group) %>%
    plot_pvalue_ital(., facet_var = "gene_symbol_vep")
  
  # Create CPG Odds Ratio plot 
  
  enr_plot <- hist_cpg_enr_all %>% 
    filter(hist_gene %in% sig_hist_cpgs_both,
           plot_group == group) %>%
    plot_enr(., facet_var = "gene_symbol_vep",
             log_scale = TRUE)
  
  # Create % patients with CPG PLP plot, and include fractions as text 
  
  perc_plot <- hist_cpg_enr_all %>%
    filter(hist_gene %in% sig_hist_cpgs_both,
           plot_group == group) %>%
    plot_perc(., facet_var = "gene_symbol_vep")
  
  # Merge plots and write to output

  merged_plot[[group]] <- ggarrange(pval_plot, enr_plot, perc_plot,
                            nrow = 1, widths = c(2.,1.25,1.65))
  
  merged_plot[[group]] <- annotate_figure(
                                          merged_plot[[group]],
                                          top = text_grob(group, face = "bold", size = 18)
                                        )

}

hist_heights <- hist_cpg_enr_all %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(cohort == "PMBB") %>%
  count(plot_group) %>%
  pull(n)

ggarrange(plotlist = merged_plot[1:length(merged_plot)],
          ncol = 1,
          heights = (1 + hist_heights/3),
          align = "hv")

ggsave(file.path(plot_dir, "sig-CPG-enrichment-PBTA-vs-control-byHist.pdf"),
        width = 9, height = 25)

write_tsv(hist_cpg_enr_all,
          file.path(results_dir, "hist-gene-plp-enr-pbta-vs-pmbb-gnomad.tsv"))
