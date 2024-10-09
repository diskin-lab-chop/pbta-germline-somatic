# Plot PBTA pathway PLP enrichment by histology
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

# Set enrichment file paths

pathway_enr_gnomad_file <- file.path(input_dir, 
                                     "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_kegg_msigdb.v2024.1.Hs.symbols_pathway_gnomad_enrichment.tsv")
pathway_enr_pmbb_file <- file.path(input_dir, 
                                   "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_kegg_msigdb.v2024.1.Hs.symbols.txt_pmbb_enrichment.tsv")

repair_enr_gnomad_file <- file.path(input_dir, 
                                    "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_dna_repair_pathway_gnomad_enrichment.tsv")
repair_enr_pmbb_file <- file.path(input_dir, 
                                  "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_dna_repair_pathway_pmbb_enrichment.tsv")

opc_hist_file <- file.path(data_dir, 
                           "histologies.tsv")

cbtn_histologies_file <- file.path(root_dir, "analyses", 
                                   "collapse-tumor-histologies", "results", 
                                   "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")
plp_all_exome_file <- file.path(root_dir, "analyses",
                                "bed-intersect", "results", 
                                "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded.tsv")

# Read in DNA repair pathway enrichment files
repair_enr_gnomad <- read_tsv(repair_enr_gnomad_file) %>%
  # rename pathways
  dplyr::mutate(pathway_name = case_when(
    grepl("Knijnenburg", pathway_name) ~ "All DNA Repair",
    pathway_name == "Mismatch Repair (MMR)" ~ "Mismatch Repair",
    pathway_name == "Others" ~ "Other DNA Repair",
    TRUE ~ pathway_name
  ))

repair_enr_pmbb <- read_tsv(repair_enr_pmbb_file) %>%
  # rename pathways
  dplyr::mutate(pathway_name = case_when(
    grepl("Knijnenburg", pathway_name) ~ "All DNA Repair",
    pathway_name == "Mismatch Repair (MMR)" ~ "Mismatch Repair",
    pathway_name == "Others" ~ "Other DNA Repair",
    TRUE ~ pathway_name
  ))

blacklisted_pathways <- read_lines(file.path(input_dir, "blacklisted-pathways.txt"))

# Read in KEGG pathway enrichment files
pathway_enr_gnomad <- read_tsv(pathway_enr_gnomad_file) %>%
  dplyr::filter(!pathway_name %in% c(blacklisted_pathways)) %>%
  dplyr::mutate(pathway_name = str_replace(pathway_name, "KEGG_MEDICUS_REFERENCE_|KEGG_MEDICUS_PATHOGEN_|KEGG_MEDICUS_VARIANT_", "MEDICUS: ")) %>%
  dplyr::mutate(pathway_name = str_replace(pathway_name, "KEGG_", "")) %>%
  dplyr::mutate(pathway_name = str_replace_all(pathway_name, "_", " ")) %>%
  arrange(pathway_name)

pathway_enr_pmbb <- read_tsv(pathway_enr_pmbb_file) %>%
  dplyr::filter(!pathway_name %in% c(blacklisted_pathways)) %>%
  dplyr::mutate(pathway_name = str_replace(pathway_name, "KEGG_MEDICUS_REFERENCE_|KEGG_MEDICUS_PATHOGEN_|KEGG_MEDICUS_VARIANT_", "MEDICUS: ")) %>%
  dplyr::mutate(pathway_name = str_replace(pathway_name, "KEGG_", "")) %>%
  dplyr::mutate(pathway_name = str_replace_all(pathway_name, "_", " ")) %>%
  arrange(pathway_name)

# Load histologies 
hist <- read_tsv(cbtn_histologies_file) %>%
  dplyr::mutate(plot_group = case_when(
    grepl("Pineoblastoma", pathology_diagnosis) ~ "Pineoblastoma",
    TRUE ~ plot_group
  ))

# calculate plot group counts for gnomad and PMBB comparisons

#obtain sample counts per plot group for gnomad comparison; here, we will filter out WXS samples
hist_cts <- hist %>%
  count(plot_group) %>%
  dplyr::rename(total_cohort_size_case = n)

# Load autogvp output for gnomad comparison 
plp <- read_tsv(plp_all_exome_file) %>%
  distinct(Kids_First_Biospecimen_ID_normal, gene_symbol_vep, .keep_all = TRUE) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_normal, plot_group))

# Get list of enrichment dfs
gnomad_list = list(pathway_enr_gnomad,
                   repair_enr_gnomad)
names(gnomad_list) <- c("KEGG_pathways", "Knijnenburg_repair_pathways")

pmbb_list = list(pathway_enr_pmbb,
                   repair_enr_pmbb)

# loop through enr dfs to calculate enrichment
for (i in 1:length(gnomad_list)){
  
  enr_gnomad <- gnomad_list[[i]]
  enr_pmbb <- pmbb_list[[i]]

  # Create list of pathways with corresponding vector of included genes
  pathway_genes <- strsplit(enr_gnomad$genes, ",")
  names(pathway_genes) <- enr_gnomad$pathway_name

  # define dfs in which to calculate enrichment
  hist_pathway_ct_df <- data.frame(plot_group = rep(unique(hist$plot_group), length(pathway_genes)),
                                   pathway_name = rep(names(pathway_genes), each = length(unique(hist$plot_group))),
                                   count_with_plp_case = 0)
  
  # get P-LP carrier counts in pathway
  for (j in 1:nrow(hist_pathway_ct_df)){
    
    group <- hist_pathway_ct_df$plot_group[j]
    pathway <- hist_pathway_ct_df$pathway_name[j]
    
    hist_pathway_ct_df$count_with_plp_case[j] <- length(unique(plp$Kids_First_Biospecimen_ID_normal[plp$plot_group %in% group & plp$gene_symbol_vep %in% pathway_genes[[pathway]]]))

  }

  # filter and append hist info to count dfs
  hist_pathway_ct_df <- hist_pathway_ct_df %>%
    dplyr::filter(count_with_plp_case > 1) %>%
    left_join(hist_cts) %>%
    dplyr::mutate(count_without_plp_case = total_cohort_size_case - count_with_plp_case)

  # Calculate P-LP carrier enrichment in pathway relative to gnomAD
  hist_enr_gnomad <- enr_gnomad %>%
    dplyr::select(pathway_name, count_with_plp_control,
                  total_cohort_size_control, count_without_plp_control) %>%
    left_join(hist_pathway_ct_df) %>%
    dplyr::filter(!is.na(plot_group)) %>%
    dplyr::mutate(p = NA_integer_,
                  OR = NA_integer_,
                  ci.int1 = NA_integer_,
                  ci.int2 = NA_integer_)
  
  hist_enr_gnomad <- calculate_enrichment(hist_enr_gnomad)

  # multiple test correction 
  hist_enr_gnomad <- hist_enr_gnomad %>%
    group_by(plot_group) %>%
    dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))
  
  # repeat for enrichment relative to PMBB
  hist_enr_pmbb <- enr_pmbb %>%
    dplyr::select(pathway_name, count_with_plp_control,
                  total_cohort_size_control, count_without_plp_control) %>%
    left_join(hist_pathway_ct_df) %>%
    dplyr::filter(!is.na(plot_group)) %>%
    dplyr::mutate(p = NA_integer_,
                  OR = NA_integer_,
                  ci.int1 = NA_integer_,
                  ci.int2 = NA_integer_)
  
  hist_enr_pmbb <- calculate_enrichment(hist_enr_pmbb)
  
  # multiple test correction
  hist_enr_pmbb <- hist_enr_pmbb %>%
    group_by(plot_group) %>%
    dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))
  
  # Create dummy PBTA df for plotting
  hist_enr_pbta <- hist_enr_pmbb %>%
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
  hist_enr_gnomad <- hist_enr_gnomad %>%
    dplyr::select(pathway_name, plot_group, count_with_plp_control, count_without_plp_control, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "gnomAD") %>%
    dplyr::rename("n_plp" = "count_with_plp_control", 
                  "n_no_plp" = "count_without_plp_control")
  
  hist_enr_pmbb <- hist_enr_pmbb %>%
    dplyr::select(pathway_name, plot_group, count_with_plp_control, count_without_plp_control, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "PMBB") %>%
    dplyr::rename("n_plp" = "count_with_plp_control", 
                  "n_no_plp" = "count_without_plp_control")
  
  # merge enrichment results and calculate percent of patients with CPG PLP
  hist_enr_all <- hist_enr_pbta %>%
    bind_rows(hist_enr_gnomad, hist_enr_pmbb) %>%
    mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
           fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
    arrange(desc(-log10(p))) %>%
    mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
           hist_pathway = glue::glue("{plot_group}: {pathway_name}"),
           hist_pathway = factor(hist_pathway, unique(hist_pathway)),) %>%
    write_tsv(file.path(results_dir, glue::glue("{names(gnomad_list)[i]}-enrichment-by-hist-PBTA-vs-control.tsv")))
  
  # pull significantly enriched pathways relative to gnomAD and PMBB cohorts, and obtain those common to both sets
  sig_hist_gnomad <- hist_enr_all %>%
    filter(cohort == "gnomAD",
           padj < 0.05,
           OR > 1) %>%
    pull(hist_pathway)
  
  sig_hist_pmbb <- hist_enr_all %>%
    filter(cohort == "PMBB",
           padj < 0.05,
           OR > 1) %>%
    pull(hist_pathway)
  
  sig_hist_both <- intersect(sig_hist_gnomad, sig_hist_pmbb)
  
  # Loop through histologies and plot significantly enriched pathways
  plot_groups <- hist_enr_all %>%
    filter(hist_pathway %in% sig_hist_both) %>%
    pull(plot_group) %>%
    unique()

  for (group in plot_groups) {
    
    n_sig <- sum(grepl(group, sig_hist_both))
    
    # Create CPG enrichment FDR plot  
    pval_plot <- hist_enr_all %>% 
      filter(plot_group == group) %>%
      plot_pvalue(., facet_var = "hist_pathway",
                  to_retain = sig_hist_both)
    
    # Create CPG Odds Ratio plot 
    enr_plot <- hist_enr_all %>% 
      filter(hist_pathway %in% sig_hist_both,
             plot_group == group) %>%
      plot_enr(., facet_var = "hist_pathway",
               log_scale = TRUE)
    
    # Create % patients with CPG PLP plot separated by source call, and include fractions as text 
    perc_plot <- hist_enr_all %>%
      filter(hist_pathway %in% sig_hist_both,
             plot_group == group) %>%
      plot_perc(., facet_var = "hist_pathway")
    
    # Merge plots and write to output
    ggarrange(pval_plot, enr_plot, perc_plot,
              nrow = 1, widths = c(2.5,1.25,1.65))
    
    ggsave(file.path(plot_dir, glue::glue("sig-{group}-{names(gnomad_list)[i]}-enrichment-PBTA-vs-control.pdf")),
           width = 10, height = 2.5 +((n_sig - 1) * 1.5))
    
  }
  
}
