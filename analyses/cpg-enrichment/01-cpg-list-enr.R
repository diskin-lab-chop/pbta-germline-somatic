# Plot PBTA CPG PLP enrichment 
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
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))
source(file.path(analysis_dir, "util", "enrichment_functions.R"))

# Set file paths

all_cpg_enr_gnomad_file <- file.path(input_dir, "pbta-merged-plp-variants-autogvp-abridged-no-wxs_cpg_pathway_gnomad_enrichment.tsv")
all_cpg_enr_pmbb_file <- file.path(input_dir, "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_cpg_pathway_pmbb_enrichment.tsv")

opc_hist_file <- file.path(data_dir, 
                           "histologies.tsv")

cbtn_histologies_file <- file.path(root_dir, "analyses", 
                                   "collapse-tumor-histologies", "results", 
                                   "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

plp_all_exome_file <- file.path(root_dir, "analyses",
                                "bed-intersect", "results", 
                                "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded.tsv")

plp_no_wxs_file <- file.path(root_dir, "analyses",
                                "bed-intersect", "results", 
                                "pbta-merged-plp-variants-autogvp-abridged-no-wxs.tsv")

cpg_file <- file.path(root_dir, "analyses", 
                      "oncokb-annotation", 
                      "input", "cpg.txt")

# Load enrichment results
all_cpg_enr_gnomad <- read_tsv(all_cpg_enr_gnomad_file) %>%
  dplyr::mutate(pathway_name = "Cancer predisposition genes")

all_cpg_enr_pmbb <- read_tsv(all_cpg_enr_pmbb_file) %>%
  dplyr::mutate(pathway_name = "Cancer predisposition genes")

# Reformat enrichment results for plotting

# create pbta enrichment results df with empty values
all_cpg_enr_pbta <- all_cpg_enr_gnomad %>%
  dplyr::select(pathway_name, count_with_plp_case, count_without_plp_case, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PBTA",
         OR = NA,
         p = NA,
         ci.int1 = NA, 
         ci.int2 = NA,
         padj = NA) %>%
  dplyr::rename("n_plp" = "count_with_plp_case", 
                "n_no_plp" = "count_without_plp_case")

# reformat gnomAD and PMBB enrichment results

cpg_enr_gnomad <- all_cpg_enr_gnomad %>%
  dplyr::select(pathway_name, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "gnomAD") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

cpg_enr_pmbb <- all_cpg_enr_pmbb %>%
  dplyr::select(pathway_name, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PMBB") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

# merge gnomAD and PMBB enrichment dfs with empty PBTA df

cpg_enr_all <- all_cpg_enr_pbta %>%
  bind_rows(cpg_enr_gnomad, cpg_enr_pmbb) %>%
  mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
         fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
  arrange(desc(-log10(p))) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")))

# Create CPG enrichment FDR plot  

pval_plot <- plot_pvalue(cpg_enr_all,
                         facet_var = "pathway_name")

# Create CPG Odds Ratio plot 

enr_plot <- plot_enr(cpg_enr_all,
                     facet_var = "pathway_name")

# Create % patients with CPG PLP plot 

perc_plot <- plot_perc(cpg_enr_all,
                       facet_var = "pathway_name")

# Merge plots and write to output

pdf(file.path(plot_dir, "all-CPG-enrichment-PBTA-vs-control.pdf"),
     width = 8, height = 2.5)

ggarrange(pval_plot, enr_plot, perc_plot,
          nrow = 1, widths = c(1.5,1.25,1.5))

dev.off()

## CPG enrichment by plot group

## gnomad enrichment

# extract all WXS IDs, to be filtered for gnomad enrichment
wxs_ids <- read_tsv(opc_hist_file) %>%
  dplyr::filter(experimental_strategy == "WXS") %>%
  pull(Kids_First_Biospecimen_ID)

# Load histologies file and obtain sample counts per plot group for gnomad comparison
hist <- read_tsv(cbtn_histologies_file)

hist_cts_gnomad <- hist %>%
  dplyr::filter(!Kids_First_Biospecimen_ID_normal %in% wxs_ids) %>%
  count(plot_group) %>%
  dplyr::rename(total_cohort_size_case = n)

# Load CPGs

cpgs <- read_lines(cpg_file)

# Load autogvp output for gnomad comparison 
plp_gnomad <- read_tsv(plp_no_wxs_file) %>%
  dplyr::filter(gene_symbol_vep %in% cpgs) %>%
  distinct(Kids_First_Biospecimen_ID_normal, .keep_all = TRUE) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_normal, plot_group))

# Obtain P-LP carrier count by plot group
hist_plp_ct_gnomad <- plp_gnomad %>%
  count(plot_group) %>%
  dplyr::rename(count_with_plp_case = n) %>%
  right_join(hist_cts_gnomad) %>%
  dplyr::mutate(count_with_plp_case = case_when(
    is.na(count_with_plp_case) ~ 0,
    TRUE ~ count_with_plp_case
  )) %>%
  dplyr::mutate(count_without_plp_case = total_cohort_size_case - count_with_plp_case) %>%
  dplyr::mutate(pathway_name = "Cancer predisposition genes")

# Calculate plot group P-LP carrier enrichment relative to gnomAD
hist_cpg_enr_gnomad <- all_cpg_enr_gnomad %>%
  dplyr::select(pathway_name, count_with_plp_control,
                total_cohort_size_control, count_without_plp_control) %>%
  left_join(hist_plp_ct_gnomad) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  dplyr::mutate(p = 0,
                OR = 0,
                ci.int1 = 0,
                ci.int2 = 0)

hist_cpg_enr_gnomad <- calculate_enrichment(hist_cpg_enr_gnomad)

# Repeat enrichment calculations relative to PMBB cohort

hist_cts_pmbb <- hist %>%
  count(plot_group) %>%
  dplyr::rename(total_cohort_size_case = n)

plp_pmbb <- read_tsv(plp_all_exome_file) %>%
  dplyr::filter(gene_symbol_vep %in% cpgs) %>%
  distinct(Kids_First_Biospecimen_ID_normal, .keep_all = TRUE) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_normal, plot_group))

# Obtain P-LP carrier count by plot group
hist_plp_ct_pmbb <- plp_pmbb %>%
  count(plot_group) %>%
  dplyr::rename(count_with_plp_case = n) %>%
  right_join(hist_cts_pmbb) %>%
  dplyr::mutate(count_with_plp_case = case_when(
    is.na(count_with_plp_case) ~ 0,
    TRUE ~ count_with_plp_case
  )) %>%
  dplyr::mutate(count_without_plp_case = total_cohort_size_case - count_with_plp_case) %>%
  dplyr::mutate(pathway_name = "Cancer predisposition genes")

hist_cpg_enr_pmbb <- all_cpg_enr_pmbb %>%
  dplyr::select(pathway_name, count_with_plp_control,
                total_cohort_size_control, count_without_plp_control) %>%
  left_join(hist_plp_ct_pmbb) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  dplyr::mutate(p = 0,
                OR = 0,
                ci.int1 = 0,
                ci.int2 = 0)

hist_cpg_enr_pmbb <- calculate_enrichment(hist_cpg_enr_pmbb)


# Create empty PBTA enrichment df for plotting
hist_cpg_enr_pbta <- hist_cpg_enr_pmbb %>%
  dplyr::select(pathway_name, plot_group, count_with_plp_case, count_without_plp_case, 
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
  dplyr::select(pathway_name, plot_group, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "gnomAD") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

hist_cpg_enr_pmbb <- hist_cpg_enr_pmbb %>%
  dplyr::select(pathway_name, plot_group, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PMBB") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

# merge PBTA, PMBB, and gnomAD dfs
hist_cpg_enr_all <- hist_cpg_enr_pbta %>%
  bind_rows(hist_cpg_enr_gnomad, hist_cpg_enr_pmbb) %>%
  mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
         fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
  arrange(desc(-log10(p))) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
         plot_group = factor(plot_group, unique(plot_group)),)

# Create CPG enrichment FDR plot faceted by plot group

hist_pval_plot <- plot_pvalue(hist_cpg_enr_all,
                              facet_var = "plot_group")
  
# Create CPG Odds Ratio plot faceted by plot group 

hist_enr_plot <- plot_enr(hist_cpg_enr_all,
                              facet_var = "plot_group")
  
# Create % patients with CPG PLP plot, faceted by plot group

hist_perc_plot <- plot_perc(hist_cpg_enr_all,
                              facet_var = "plot_group")
  
  
# Merge plots and write to output
  
ggarrange(hist_pval_plot, hist_enr_plot, hist_perc_plot,
          nrow = 1, widths = c(1.7,1.25,1.65))

ggsave(file.path(plot_dir, glue::glue("hist-all-CPG-enrichment-PBTA-vs-control.pdf")),
       width = 8, height = 22)
  