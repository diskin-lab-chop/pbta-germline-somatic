# Plot PBTA CPG PLP enrichment by plot group
# Ryan Corbett
# F2024

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
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))
source(file.path(analysis_dir, "util", "enrichment_functions.R"))

# Set file file paths

cpg_enr_gnomad_file <- file.path(input_dir, 
                                 "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_gnomAD_enrichment.tsv")
cpg_enr_pmbb_file <- file.path(input_dir, 
                               "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_PMBB_enrichment.tsv")

cbtn_histologies_file <- file.path(root_dir, "analyses", 
                                   "collapse-tumor-histologies", "results", 
                                   "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")
plp_file <- file.path(data_dir, 
                      "pbta-merged-plp-variants-autogvp-abridged.tsv")

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

# get plot groups counts
hist_cts <- hist %>%
  count(plot_group) %>%
  dplyr::rename(total_cohort_size_case = n)

# Read in autogvp results and filter out intronic NF1 plp variants
plp <- read_tsv(plp_file) %>%
  distinct(Kids_First_Biospecimen_ID_normal, gene_symbol_vep, .keep_all = TRUE) %>%
  dplyr::filter(gene_symbol_vep %in% cpgs,
                !grepl("Koczkowska", autogvp_call_reason)) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_normal, plot_group))

# Obtain P-LP carrier counts by plot group
hist_plp_ct <- plp %>%
  count(gene_symbol_vep, plot_group) %>%
  dplyr::filter(n > 1) %>%
  dplyr::rename(count_with_plp_case = n) %>%
  left_join(hist_cts) %>%
  dplyr::mutate(count_without_plp_case = total_cohort_size_case - count_with_plp_case)
  
# Calculate plot-group level P-LP carrier enrichment relative to gnomAD
hist_cpg_enr_gnomad <- cpg_enr_gnomad %>%
  dplyr::filter(!is.na(gene_symbol_vep)) %>%
  dplyr::select(gene_symbol_vep, count_with_plp_control,
                total_cohort_size_control, count_without_plp_control) %>%
  left_join(hist_plp_ct) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  dplyr::mutate(p = 0,
                OR = 0,
                ci.int1 = 0,
                ci.int2 = 0)

for (i in 1:nrow(hist_cpg_enr_gnomad)){
  
  fisher_test <- fisher.test(matrix(c(hist_cpg_enr_gnomad$count_with_plp_case[i],
                                      hist_cpg_enr_gnomad$count_without_plp_case[i],
                                      hist_cpg_enr_gnomad$count_with_plp_control[i],
                                      hist_cpg_enr_gnomad$count_without_plp_control[i]),
                                    2, 2))
  
  hist_cpg_enr_gnomad$OR[i] = fisher_test$estimate
  hist_cpg_enr_gnomad$p[i] = fisher_test$p.value
  hist_cpg_enr_gnomad$ci.int1[i] = fisher_test$conf.int[1]
  hist_cpg_enr_gnomad$ci.int2[i] = fisher_test$conf.int[2]
  
}

# multiple test correction 

hist_cpg_enr_gnomad <- hist_cpg_enr_gnomad %>%
  group_by(plot_group) %>%
  dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))

# Repeat enrichment calculations for PMBB controls

hist_cpg_enr_pmbb <- cpg_enr_pmbb %>%
  dplyr::filter(!is.na(gene_symbol_vep)) %>%
  dplyr::select(gene_symbol_vep, count_with_plp_control,
                total_cohort_size_control, count_without_plp_control) %>%
  left_join(hist_plp_ct) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  dplyr::mutate(p = 0,
                OR = 0,
                ci.int1 = 0,
                ci.int2 = 0)

for (i in 1:nrow(hist_cpg_enr_gnomad)){
  
  fisher_test <- fisher.test(matrix(c(hist_cpg_enr_pmbb$count_with_plp_case[i],
                                      hist_cpg_enr_pmbb$count_without_plp_case[i],
                                      hist_cpg_enr_pmbb$count_with_plp_control[i],
                                      hist_cpg_enr_pmbb$count_without_plp_control[i]),
                                    2, 2))
  
  hist_cpg_enr_pmbb$OR[i] = fisher_test$estimate
  hist_cpg_enr_pmbb$p[i] = fisher_test$p.value
  hist_cpg_enr_pmbb$ci.int1[i] = fisher_test$conf.int[1]
  hist_cpg_enr_pmbb$ci.int2[i] = fisher_test$conf.int[2]
  
}

# mulitiple test correction

hist_cpg_enr_pmbb <- hist_cpg_enr_pmbb %>%
  group_by(plot_group) %>%
  dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))

# Create 'dummy' data frame with CBTN enrichment values to use as reference in plots

hist_cpg_enr_pbta <- hist_cpg_enr_gnomad %>%
  dplyr::select(gene_symbol_vep, plot_group, count_with_plp_case, count_without_plp_case, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PBTA",
         OR = NA,
         p = NA,
         ci.int1 = NA, 
         ci.int2 = NA,
         padj = NA) %>%
  dplyr::rename("n_plp" = "count_with_plp_case", 
                "n_no_plp" = "count_without_plp_case")

# subset gnomAD and PMBB enrichment results for merging

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
         hist_gene = glue::glue("{gene_symbol_vep}:{plot_group}"),
         hist_gene = factor(hist_gene, unique(hist_gene)),)

# pull significnalty enriched CPGs relative to gnomAD and PMBB cohorts, and obtain those common to both sets

sig_hist_cpgs_gnomad <- hist_cpg_enr_all %>%
  filter(cohort == "gnomAD" & padj < 0.05 & padj > 0) %>%
  pull(hist_gene)

sig_hist_cpgs_pmbb <- hist_cpg_enr_all %>%
  filter(cohort == "PMBB" & padj < 0.05 & padj > 0) %>%
  pull(hist_gene)

sig_hist_cpgs_both <- intersect(sig_hist_cpgs_gnomad, sig_hist_cpgs_pmbb)

# Loop through plot groups to plot sig results
plot_groups <- hist_cpg_enr_all %>%
  filter(hist_gene %in% sig_hist_cpgs_both) %>%
  pull(plot_group) %>%
  unique()

for (group in plot_groups) {
  
  # get number of significantly enriched genes in plot group to determine plot height
  n_sig <- sum(grepl(group, sig_hist_cpgs_both))
  
  # Create CPG enrichment FDR plot  
  
  pval_plot <- hist_cpg_enr_all %>%
    dplyr::filter(plot_group == group) %>%
    plot_pvalue(., facet_var = "hist_gene",
                to_retain = sig_hist_cpgs_both)
  
  # Create CPG Odds Ratio plot 
  
  enr_plot <- hist_cpg_enr_all %>% 
    filter(hist_gene %in% sig_hist_cpgs_both,
           plot_group == group) %>%
    plot_enr(., facet_var = "hist_gene",
             log_scale = TRUE)
  
  # Create % patients with CPG PLP plot, and include fractions as text 
  
  perc_plot <- hist_cpg_enr_all %>%
    filter(hist_gene %in% sig_hist_cpgs_both,
           plot_group == group) %>%
    plot_perc(., facet_var = "hist_gene")
  
  # Merge plots and write to output

  ggarrange(pval_plot, enr_plot, perc_plot,
            nrow = 1, widths = c(1.75,1.5,1.5))
  
  ggsave(file.path(plot_dir, glue::glue("sig-{group}-CPG-enrichment-PBTA-vs-control.tiff")),
         width = 9, height = 2.5 +((n_sig - 1) * 1.5), units = "in")

}
