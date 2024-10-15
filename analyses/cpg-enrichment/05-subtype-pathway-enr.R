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

# Set enrichment file paths

cpg_enr_gnomad_file <- file.path(input_dir, "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_gnomAD_enrichment.tsv")
cpg_enr_pmbb_file <- file.path(input_dir, "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_PMBB_enrichment.tsv")

pathway_enr_gnomad_file <- file.path(input_dir, "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_kegg_gnomAD_enrichment.tsv")
pathway_enr_pmbb_file <- file.path(input_dir, "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_kegg_PMBB_enrichment.tsv")

repair_enr_gnomad_file <- file.path(input_dir, "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_Knijnenburg_gnomAD_enrichment.tsv")
repair_enr_pmbb_file <- file.path(input_dir, "pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_Knijnenburg_PMBB_enrichment.tsv")

#cbtn_histologies_file <- file.path(root_dir,"analyses", "collapse-tumor-histologies", "results", "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")
cbtn_histologies_file <- file.path(root_dir,"analyses", "survival", "results", "germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv")

plp_file <- file.path(data_dir, "pbta-merged-plp-variants-autogvp-abridged.tsv")

# Read in hist, plp, and enrichment by CPG files



repair_enr_gnomad <- read_tsv(repair_enr_gnomad_file) %>%
  dplyr::mutate(pathway_name = case_when(
    pathway_name == "Mismatch Repair (MMR)" ~ "Mismatch Repair",
    pathway_name == "Others" ~ "Other DNA Repair",
    TRUE ~ pathway_name
  ))

repair_enr_pmbb <- read_tsv(repair_enr_pmbb_file) %>%
  dplyr::mutate(pathway_name = case_when(
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


hist <- read_tsv(cbtn_histologies_file)

hist_cts <- hist %>%
  count(mol_sub_group) %>%
  dplyr::rename(total_cohort_size_case = n)

plp <- read_tsv(plp_file) %>%
  distinct(Kids_First_Biospecimen_ID_normal, gene_symbol_vep, .keep_all = TRUE) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_normal, mol_sub_group))

gnomad_list = list(pathway_enr_gnomad,
                   repair_enr_gnomad)
names(gnomad_list) <- c("KEGG_pathways", "Knijnenburg_repair_pathways")

pmbb_list = list(pathway_enr_pmbb,
                   repair_enr_pmbb)

for (i in 1:length(gnomad_list)){
  
  enr_gnomad <- gnomad_list[[i]]
  enr_pmbb <- pmbb_list[[i]]

  pathway_genes <- strsplit(enr_gnomad$genes, ",")
  names(pathway_genes) <- enr_gnomad$pathway_name

  hist_pathway_ct_df <- data.frame(mol_sub_group = rep(unique(hist$mol_sub_group), length(pathway_genes)),
                                   pathway_name = rep(names(pathway_genes), each = length(unique(hist$mol_sub_group))),
                                   count_with_plp_case = 0)

  for (j in 1:nrow(hist_pathway_ct_df)){
    
    group <- hist_pathway_ct_df$mol_sub_group[j]
    pathway <- hist_pathway_ct_df$pathway_name[j]
    
    hist_pathway_ct_df$count_with_plp_case[j] <- length(unique(plp$Kids_First_Biospecimen_ID_normal[plp$mol_sub_group %in% group & plp$gene_symbol_vep %in% pathway_genes[[pathway]]]))
    
  }

  hist_pathway_ct_df <- hist_pathway_ct_df %>%
    dplyr::filter(count_with_plp_case > 1) %>%
    left_join(hist_cts) %>%
    dplyr::mutate(count_without_plp_case = total_cohort_size_case - count_with_plp_case)

  hist_enr_gnomad <- enr_gnomad %>%
    dplyr::select(pathway_name, count_with_plp_control,
                  total_cohort_size_control, count_without_plp_control) %>%
    left_join(hist_pathway_ct_df) %>%
    dplyr::filter(!is.na(mol_sub_group)) %>%
    dplyr::mutate(p = 0,
                  OR = 0,
                  ci.int1 = 0,
                  ci.int2 = 0)

  for (j in 1:nrow(hist_enr_gnomad)){
    
    fisher_test <- fisher.test(matrix(c(hist_enr_gnomad$count_with_plp_case[j],
                                        hist_enr_gnomad$count_without_plp_case[j],
                                        hist_enr_gnomad$count_with_plp_control[j],
                                        hist_enr_gnomad$count_without_plp_control[j]),
                                      2, 2))
    
    hist_enr_gnomad$OR[j] = fisher_test$estimate
    hist_enr_gnomad$p[j] = fisher_test$p.value
    hist_enr_gnomad$ci.int1[j] = fisher_test$conf.int[1]
    hist_enr_gnomad$ci.int2[j] = fisher_test$conf.int[2]
    
  }
  
  hist_enr_gnomad <- hist_enr_gnomad %>%
    group_by(mol_sub_group) %>%
    dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))
  
  
  hist_enr_pmbb <- enr_pmbb %>%
    dplyr::select(pathway_name, count_with_plp_control,
                  total_cohort_size_control, count_without_plp_control) %>%
    left_join(hist_pathway_ct_df) %>%
    dplyr::filter(!is.na(mol_sub_group)) %>%
    dplyr::mutate(p = 0,
                  OR = 0,
                  ci.int1 = 0,
                  ci.int2 = 0)
  
  for (j in 1:nrow(hist_enr_pmbb)){
    
    fisher_test <- fisher.test(matrix(c(hist_enr_pmbb$count_with_plp_case[j],
                                        hist_enr_pmbb$count_without_plp_case[j],
                                        hist_enr_pmbb$count_with_plp_control[j],
                                        hist_enr_pmbb$count_without_plp_control[j]),
                                      2, 2))
    
    hist_enr_pmbb$OR[j] = fisher_test$estimate
    hist_enr_pmbb$p[j] = fisher_test$p.value
    hist_enr_pmbb$ci.int1[j] = fisher_test$conf.int[1]
    hist_enr_pmbb$ci.int2[j] = fisher_test$conf.int[2]
    
  }
  
  hist_enr_pmbb <- hist_enr_pmbb %>%
    group_by(mol_sub_group) %>%
    dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))
  
  
  hist_enr_pbta <- hist_enr_gnomad %>%
    dplyr::select(pathway_name, mol_sub_group, count_with_plp_case, count_without_plp_case, 
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
    dplyr::select(pathway_name, mol_sub_group, count_with_plp_control, count_without_plp_control, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "gnomAD") %>%
    dplyr::rename("n_plp" = "count_with_plp_control", 
                  "n_no_plp" = "count_without_plp_control")
  
  hist_enr_pmbb <- hist_enr_pmbb %>%
    dplyr::select(pathway_name, mol_sub_group, count_with_plp_control, count_without_plp_control, 
                  OR, p, ci.int1, ci.int2, padj) %>%
    mutate(cohort = "PMBB") %>%
    dplyr::rename("n_plp" = "count_with_plp_control", 
                  "n_no_plp" = "count_without_plp_control")
  
  hist_enr_all <- hist_enr_pbta %>%
    bind_rows(hist_enr_gnomad, hist_enr_pmbb) %>%
    mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
           fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
    arrange(desc(-log10(p))) %>%
    mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
           hist_pathway = glue::glue("{pathway_name}:{mol_sub_group}"),
           hist_pathway = factor(hist_pathway, unique(hist_pathway)),) %>%
    write_tsv(file.path(results_dir, glue::glue("{names(gnomad_list)[i]}-enrichment-by-hist-PBTA-vs-control.tsv")))
  
  # pull significnalty enriched CPGs relative to gnomAD and PMBB cohorts, and obtain those common to both sets
  
  sig_hist_gnomad <- hist_enr_all %>%
    filter(cohort == "gnomAD",
           padj < 0.05,
           OR > 1,
           !grepl("cancer|carcinoma|leukemia|Melanoma|carcinogenesis|Glioma", hist_pathway)) %>%
    pull(hist_pathway)
  
  sig_hist_pmbb <- hist_enr_all %>%
    filter(cohort == "PMBB",
           padj < 0.05,
           OR > 1,
           !grepl("cancer|carcinoma|leukemia|Melanoma|carcinogenesis|Glioma", hist_pathway)) %>%
    pull(hist_pathway)
  
  sig_hist_both <- intersect(sig_hist_gnomad, sig_hist_pmbb)
  
  
  mol_sub_groups <- hist_enr_all %>%
    filter(hist_pathway %in% sig_hist_both) %>%
    pull(mol_sub_group) %>%
    unique()

  for (group in mol_sub_groups) {
    
    n_sig <- sum(grepl(group, sig_hist_both))
    pval_x_max <- hist_enr_all %>%
      dplyr::filter(mol_sub_group == group,
                    hist_pathway %in% sig_hist_both,
                    !is.na(p)) %>%
      pull(p) %>%
      min()
    
    perc_x_max <- hist_enr_all %>%
      dplyr::filter(mol_sub_group == group,
                    hist_pathway %in% sig_hist_both) %>%
      pull(perc_plp) %>%
      max()
    
    # Create CPG enrichment FDR plot  
    
    pval_plot <- hist_enr_all %>% 
      filter(hist_pathway %in% sig_hist_both,
             mol_sub_group == group) %>%
      mutate(p = -log10(p)) %>%
      ggplot(aes(x = p, y = factor(cohort))) +
      geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
      labs(x = "-log10(p)", y = "") + 
      xlim(0, -log10(pval_x_max)*1.2) +
      geom_vline(xintercept = -log10(0.05/(sum(hist_enr_all$mol_sub_group == group)/3)), linetype = "dashed") + 
      facet_wrap(~hist_pathway, nrow = length(sig_hist_both), scale = "fixed",
                 labeller = labeller(hist_pathway = label_wrap_gen(18))) + 
      theme_Publication() +
      theme(plot.margin = unit(c(2,1,1,0), "lines"))
    
    
    # Create CPG Odds Ratio plot 
    
    enr_plot <- hist_enr_all %>% 
      filter(hist_pathway %in% sig_hist_both,
             mol_sub_group == group) %>%
      ggplot(aes(x = log10(OR), y = factor(cohort))) +
      geom_point(size = 3, color = "#00A087FF",
                 show.legend = FALSE) + 
      geom_errorbar(aes(xmin = log10(ci.int1), xmax = log10(ci.int2)), width = 0.2, 
                    show.legend = FALSE, color = "#00A087FF") +
      labs(x = "log10-OR (95% CI)", y = NULL) + 
      xlim(-1, NA) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_y_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                                "PMBB" = "")) +
      facet_wrap(~hist_pathway, nrow = length(sig_hist_both), scale = "fixed",
                 labeller = labeller(hist_pathway = label_wrap_gen(18))) +  expand_limits(y=0) +
      theme_Publication() +
      theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"))
    
    
    
    # Create % patients with CPG PLP plot separated by source call, and include fractions as text 
    
    perc_plot <- hist_enr_all %>%
      filter(hist_pathway %in% sig_hist_both,
             mol_sub_group == group) %>%
      ggplot(aes(x = perc_plp, y = factor(cohort), label = fraction)) +
      geom_bar(stat = "identity", color = "black",
               fill = "#00A087FF",
               show.legend = TRUE) + 
      geom_text(x = perc_x_max * 1.05, hjust = 0, size = 4, fontface = 2) +
      labs(x = "% Cohort P/LP", y = NULL, fill = NULL) + 
      xlim(0, perc_x_max) +
      scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                                "PMBB" = NULL)) +
      guides(fill = guide_legend(nrow = 1)) +
      facet_wrap(~hist_pathway, nrow = length(sig_hist_both), scale = "fixed",
                 labeller = labeller(hist_pathway = label_wrap_gen(18))) +
      expand_limits(x=70) +
      coord_cartesian(clip = 'off') +
      theme_Publication() +
      theme(plot.margin = unit(c(2,4.5,1,1), "lines"))
    
    
    # Merge plots and write to output
    
    ggarrange(pval_plot, enr_plot, perc_plot,
              nrow = 1, widths = c(1.7,1.25,1.65))
    
    ggsave(file.path(plot_dir, glue::glue("sig-{group}-{names(gnomad_list)[i]}-enrichment-PBTA-vs-control.tiff")),
           width = 9, height = 2.5 +((n_sig - 1) * 2), units = "in")
    
  }
  
}

level_order = hist_enr_all %>%
  filter(pathway_name == "DNA repair",
         cohort == "PMBB") %>%
  arrange(desc(p)) %>%
  pull(mol_sub_group)

pval_plot <- hist_enr_all %>%
  filter(pathway_name == "DNA repair",
         cohort == "PMBB") %>%
  dplyr::mutate(mol_sub_group = as.character(mol_sub_group)) %>%
  dplyr::mutate(p = -log10(p),
         mol_sub_group = factor(mol_sub_group, levels = level_order)) %>%

  ggplot(aes(x = p, y = factor(mol_sub_group))) +
  geom_point(size = 3,
             color = "#00A087FF",
             show.legend = FALSE) + 
  labs(x = "-log10(p) DNA Repair P/LP enrichment", y = "") + 
  xlim(0, 5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") + 
  theme_Publication() +
  theme(plot.margin = unit(c(2,1,1,0), "lines"))

pval_plot

ggsave(file.path(plot_dir, "dna-repair-subgroup-enr-pmbb.pdf"),
       width = 6, height = 6)
