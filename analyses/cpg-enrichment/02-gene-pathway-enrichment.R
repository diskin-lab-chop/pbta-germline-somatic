# Plot PBTA CPG PLP enrichment 
# Ryan Corbett
# February 2023


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

# Set enrichment file paths

cpg_enr_gnomad_file <- file.path(input_dir, 
                                 "pbta-merged-plp-variants-autogvp-abridged-no-wxs_gene_gnomad_enrichment.tsv")
cpg_enr_pmbb_file <- file.path(input_dir,
                               "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_gene_pmbb_enrichment.tsv")

pathway_enr_gnomad_file <- file.path(input_dir, 
                                     "pbta-merged-plp-variants-autogvp-abridged-no-wxs_kegg_pathway_gnomad_enrichment.tsv")
pathway_enr_pmbb_file <- file.path(input_dir, 
                                   "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_kegg_pathway_pmbb_enrichment.tsv")

repair_enr_gnomad_file <- file.path(input_dir, 
                                    "pbta-merged-plp-variants-autogvp-abridged-no-wxs_dna_repair_pathway_gnomad_enrichment.tsv")
repair_enr_pmbb_file <- file.path(input_dir, 
                                  "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_dna_repair_pathway_pmbb_enrichment.tsv")

cbtn_histologies_file <- file.path(root_dir, "analyses", 
                                   "collapse-tumor-histologies", "results", 
                                   "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

plp_all_exome_file <- file.path(root_dir, "analyses",
                                "bed-intersect", "results", 
                                "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded.tsv")

plp_no_wxs_file <- file.path(root_dir, "analyses",
                             "bed-intersect", "results", 
                             "pbta-merged-plp-variants-autogvp-abridged-no-wxs.tsv")

# Read in hist, plp, and enrichment by CPG files

cpg_enr_gnomad <- read_tsv(cpg_enr_gnomad_file)
cpg_enr_pmbb <- read_tsv(cpg_enr_pmbb_file)

hist <- read_tsv(cbtn_histologies_file)

# Create 'dummy' data frame with PBTA enrichment values to use as reference in plots

cpg_enr_pbta <- cpg_enr_pmbb %>%
  dplyr::select(gene_symbol_vep, count_with_plp_case, count_without_plp_case, 
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

cpg_enr_gnomad <- cpg_enr_gnomad %>%
  dplyr::select(gene_symbol_vep, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "gnomAD") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

cpg_enr_pmbb <- cpg_enr_pmbb %>%
  dplyr::select(gene_symbol_vep, count_with_plp_control, count_without_plp_control, 
                OR, p, ci.int1, ci.int2, padj) %>%
  mutate(cohort = "PMBB") %>%
  dplyr::rename("n_plp" = "count_with_plp_control", 
                "n_no_plp" = "count_without_plp_control")

# merge enrichment results and calculate percent of patients with CPG PLP

cpg_enr_all <- cpg_enr_pbta %>%
  bind_rows(cpg_enr_gnomad, cpg_enr_pmbb) %>%
  mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
         fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
  arrange(desc(-log10(p))) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
         gene_symbol_vep = factor(gene_symbol_vep, unique(gene_symbol_vep)))

# pull significnalty enriched CPGs relative to gnomAD and PMBB cohorts, and obtain those common to both sets

sig_cpgs_gnomad <- cpg_enr_all %>%
  filter(cohort == "gnomAD" & padj < 0.05 & padj > 0) %>%
  pull(gene_symbol_vep)

sig_cpgs_pmbb <- cpg_enr_all %>%
  filter(cohort == "PMBB" & padj < 0.05 & padj > 0) %>%
  pull(gene_symbol_vep)

sig_cpgs_both <- intersect(sig_cpgs_gnomad, sig_cpgs_pmbb)

# Create enrichment FDR plot of sig enriched CPGS 

pval_plot <- cpg_enr_all %>% 
  plot_pvalue(., facet_var = "gene_symbol_vep",
              to_retain = sig_cpgs_both)

# Create Odds Ratio plot of sig enriched CPGs

enr_plot <- cpg_enr_all %>% 
  filter(gene_symbol_vep %in% sig_cpgs_both) %>%
  plot_enr(., facet_var = "gene_symbol_vep",
                     log_scale = TRUE)

# Create % patients with CPG PLP plot separated by source call, and include fractions as text 

perc_plot <- cpg_enr_all %>% 
  filter(gene_symbol_vep %in% sig_cpgs_both) %>%
  plot_perc(., facet_var = "gene_symbol_vep")

# Merge plots and write to output

pdf(file.path(plot_dir, "sig-CPG-enrichment-PBTA-vs-control.pdf"),
     width = 9, height = 6)

ggarrange(pval_plot, enr_plot, perc_plot,
          nrow = 1, widths = c(2,1.25,1.75))

dev.off()


## KEGG and DNA Repair pathway enrichment

# Wrangle data

repair_enr_gnomad <- read_tsv(repair_enr_gnomad_file) %>%
  dplyr::mutate(pathway_name = case_when(
    grepl("Knijnenburg", pathway_name) ~ "All DNA Repair",
    pathway_name == "Mismatch Repair (MMR)" ~ "Mismatch Repair",
    pathway_name == "Others" ~ "Other DNA Repair",
    TRUE ~ pathway_name
  ))

repair_enr_pmbb <- read_tsv(repair_enr_pmbb_file) %>%
  dplyr::mutate(pathway_name = case_when(
    grepl("Knijnenburg", pathway_name) ~ "All DNA Repair",
    pathway_name == "Mismatch Repair (MMR)" ~ "Mismatch Repair",
    pathway_name == "Others" ~ "Other DNA Repair",
    TRUE ~ pathway_name
  ))

pathway_enr_gnomad <- read_tsv(pathway_enr_gnomad_file) %>%
  dplyr::filter(!grepl("multiple species|other", pathway_name)) %>%
  dplyr::mutate(pathway_name = str_split(pathway_name, " -", simplify = TRUE)[,1])

pathway_enr_pmbb <- read_tsv(pathway_enr_pmbb_file) %>%
  dplyr::filter(!grepl("multiple species|other", pathway_name)) %>%
  dplyr::mutate(pathway_name = str_split(pathway_name, " -", simplify = TRUE)[,1])

# add KEGG pathway and DNA repair pathway results to df lists
gnomad_list = list(pathway_enr_gnomad,
                   repair_enr_gnomad)
names(gnomad_list) <- c("KEGG_pathways", "Knijnenburg_repair_pathways")

pmbb_list = list(pathway_enr_pmbb,
                 repair_enr_pmbb)

# Loop through enr df list
for (i in 1:length(gnomad_list)){
  
  enr_gnomad <- gnomad_list[[i]]
  enr_pmbb <- pmbb_list[[i]]

  # create dummy df of PBTA enrichment for plotting
  pathway_enr_pbta <- enr_pmbb %>%
    dplyr::select(pathway_name, count_with_plp_case, count_without_plp_case, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "PBTA",
           OR = NA_integer_,
           p = NA_integer_,
           ci.int1 = NA_integer_, 
           ci.int2 = NA_integer_,
           padj = NA_integer_) %>%
    dplyr::rename("n_plp" = "count_with_plp_case", 
                  "n_no_plp" = "count_without_plp_case")

  # subset gnomAD and PMBB enrichment results for merging
  pathway_enr_gnomad <- enr_gnomad %>%
    dplyr::select(pathway_name, count_with_plp_control, count_without_plp_control, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "gnomAD") %>%
    dplyr::rename("n_plp" = "count_with_plp_control", 
                  "n_no_plp" = "count_without_plp_control")

  pathway_enr_pmbb <- enr_pmbb %>%
    dplyr::select(pathway_name, count_with_plp_control, count_without_plp_control, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "PMBB") %>%
    dplyr::rename("n_plp" = "count_with_plp_control", 
                  "n_no_plp" = "count_without_plp_control")

  # merge enrichment results and calculate percent of patients with CPG PLP
  pathway_enr_all <- pathway_enr_pbta %>%
    bind_rows(pathway_enr_gnomad, pathway_enr_pmbb) %>%
    mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
           fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
    arrange(desc(-log10(p))) %>%
    mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
           pathway_name = factor(pathway_name, unique(pathway_name)))

  # pull significantly enriched CPGs relative to gnomAD and PMBB cohorts, and obtain those common to both sets
  sig_pathways_gnomad <- pathway_enr_all %>%
    filter(cohort == "gnomAD" & padj < 0.05 & padj > 0) %>%
    pull(pathway_name)
  
  sig_pathways_pmbb <- pathway_enr_all %>%
    filter(cohort == "PMBB" & padj < 0.05 & padj > 0) %>%
    pull(pathway_name)

  sig_pathways_both <- intersect(sig_pathways_gnomad, sig_pathways_pmbb)
  
  # filter out cancer pathways 
  sig_pathways_both <- sig_pathways_both[!grepl("cancer|carcinoma|leukemia|Melanoma|carcinogenesis|Glioma", sig_pathways_both)]

  # Create enrichment FDR plot of sig enriched pathways
  pval_plot <- pathway_enr_all %>%
    plot_pvalue(., facet_var = "pathway_name",
                to_retain = sig_pathways_both)
  
  # Create Odds Ratio plot of sig enriched pathways
  enr_plot <- pathway_enr_all %>%
    filter(pathway_name %in% sig_pathways_both) %>%
    plot_enr(., facet_var = "pathway_name",
             log_scale = FALSE)
  
  # Create % patients with PLP plot of sig enriched pathways
  perc_plot <- pathway_enr_all %>%
    filter(pathway_name %in% sig_pathways_both) %>%
    plot_perc(., facet_var = "pathway_name")

  # Merge plots and write to output
  ggarrange(pval_plot, enr_plot, perc_plot,
            nrow = 1, widths = c(2.5,1.25,1.65))
  
  ggsave(file.path(plot_dir, glue::glue("sig-{names(gnomad_list)[i]}-enrichment-PBTA-vs-control.pdf")),
       width = 10, height = 2.5 + ((length(sig_pathways_both)-1) * 1.15))
  
  write_tsv(pathway_enr_all,
            file.path(results_dir, glue::glue("{names(gnomad_list)[i]}-enrichment-pbta-vs-pmbb-gnomad.tsv")))
  
}