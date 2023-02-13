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
input_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

# Set enrichment file paths

cpg_enr_gnomad_file <- file.path(input_dir, "CBTN_CPG_PLP_enrichment_gnomAD.tsv")
cpg_enr_pmbb_file <- file.path(input_dir, "CBTN_CPG_PLP_enrichment_PMBB.tsv")
cpg_hist_enr_gnomad_file <- file.path(input_dir, "CBTN_CPG_PLP_enrichment_by_plotGroup_gnomAD.tsv")
cpg_hist_enr_pmbb_file <- file.path(input_dir, "CBTN_CPG_PLP_enrichment_by_plotGroup_PMBB.tsv")

cbtn_histologies_file <- file.path(root_dir,"analyses", "collapse-tumor-histologies", "results", "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")
plp_file <- file.path(data_dir, "v3", "pbta_germline_plp_calls.tsv")

# Read in hist, plp, and enrichment by CPG files

cpg_enr_gnomad <- read_tsv(cpg_enr_gnomad_file)
cpg_enr_pmbb <- read_tsv(cpg_enr_pmbb_file)

hist <- read_tsv(cbtn_histologies_file)

plp <- read_tsv(plp_file) %>%
  filter(Hugo_Symbol %in% cpg_enr_gnomad$Hugo_symbol & Kids_First_Biospecimen_ID %in% hist$Kids_First_Biospecimen_ID_normal)

# Create 'dummy' data frame with CBTN enrichment values to use as reference in plots

cpg_enr_pbta <- cpg_enr_gnomad %>%
  select(Hugo_symbol, n, n_without, oddsRatio, pvalue, CI_lower, CI_upper, adjusted_p) %>%
  mutate(cohort = "PBTA",
         oddsRatio = NA,
         pvalue = NA,
         CI_lower = NA, 
         CI_upper = NA,
         adjusted_p = NA) %>%
  rename("n_plp" = n, "n_no_plp" = n_without)

# subset gnomAD and PMBB enrichment results for merging

cpg_enr_gnomad <- cpg_enr_gnomad %>%
  select(Hugo_symbol, count_with_plp_gnomAD, count_without_plp_gnomAD, 
         oddsRatio, pvalue, CI_lower, CI_upper, adjusted_p) %>%
  mutate(cohort = "gnomAD") %>%
  rename("n_plp" = count_with_plp_gnomAD, 
         "n_no_plp" = count_without_plp_gnomAD)

cpg_enr_pmbb <- cpg_enr_pmbb %>%
  select(Hugo_symbol, count_with_plp_PMBB_noTumor, 
         count_without_plp_PMBB_noTumor, 
         oddsRatio, pvalue, CI_lower, CI_upper, adjusted_p) %>%
  mutate(cohort = "PMBB") %>%
  rename("n_plp" = count_with_plp_PMBB_noTumor, 
         "n_no_plp" = count_without_plp_PMBB_noTumor)

# merge enrichment results and calculate percent of patients with CPG PLP

cpg_enr_all <- cpg_enr_pbta %>%
  bind_rows(cpg_enr_gnomad, cpg_enr_pmbb) %>%
  mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
         fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/")) %>%
  arrange(desc(-log10(pvalue))) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
         Hugo_symbol = factor(Hugo_symbol, unique(Hugo_symbol)))

# pull significnalty enriched CPGs relative to gnomAD and PMBB cohorts, and obtain those common to both sets

sig_cpgs_gnomad <- cpg_enr_all %>%
  filter(cohort == "gnomAD" & adjusted_p < 0.05 & adjusted_p > 0) %>%
  pull(Hugo_symbol)

sig_cpgs_pmbb <- cpg_enr_all %>%
  filter(cohort == "PMBB" & adjusted_p < 0.05 & adjusted_p > 0) %>%
  pull(Hugo_symbol)

sig_cpgs_both <- intersect(sig_cpgs_gnomad, sig_cpgs_pmbb)


# Create table of PBTA CPG PLP counts, separated by those called by ClinVar and Intervar

pbta_source_cts <- plp %>%
  
  # Classify call source as clinVar or interVar 
  mutate(final_call_source = case_when(grepl("ClinVar", Reasoning_for_call) ~ "ClinVar",
                                       grepl("InterVar", Reasoning_for_call) ~ "InterVar")) %>%
  
  # Get CPG PLP counts by gene and call source
  count(Hugo_Symbol, final_call_source) %>%
  rename("n_plp" = n,
         "Hugo_symbol" = Hugo_Symbol) %>%
  mutate(n_no_plp = 838-n_plp,
           perc_plp = n_plp/838*100,
         cohort = "PBTA")

# Merge this with control cohort data

all_source_cts <- cpg_enr_all %>%
  filter(cohort != "PBTA") %>%
  mutate(final_call_source = NA) %>%
  dplyr::select(Hugo_symbol, cohort, final_call_source, n_plp, n_no_plp, perc_plp) %>%
  bind_rows(pbta_source_cts) %>%
  left_join(cpg_enr_all[,c("Hugo_symbol", "cohort", "fraction")], by = c("Hugo_symbol", "cohort")) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
         Hugo_symbol = factor(Hugo_symbol, unique(cpg_enr_all$Hugo_symbol)))


# Create CPG enrichment FDR plot  

pval_plot <- cpg_enr_all %>% 
  filter(Hugo_symbol %in% sig_cpgs_both) %>%
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = pvalue, y = factor(cohort))) +
  geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
  labs(x = "-log10(p)", y = "") + 
  geom_vline(xintercept = -log10(0.05/nrow(cpg_enr_all)), linetype = "dashed") + 
  facet_wrap(~Hugo_symbol, nrow = 6, scale = "fixed") +
  coord_cartesian(xlim = c(0,24)) + 
  theme_Publication() +
  theme(plot.margin = unit(c(2,1,1,0), "lines"))



# Create CPG Odds Ratio plot 

enr_plot <- cpg_enr_all %>% 
  filter(Hugo_symbol %in% sig_cpgs_both) %>%
  ggplot(aes(x = factor(cohort), y = oddsRatio)) +
  geom_point(size = 3, color = "#00A087FF",
             show.legend = FALSE) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, 
                show.legend = FALSE, color = "#00A087FF") +
  labs(y = "Odds Ratio (95% CI)", x = NULL) + 
  scale_x_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                            "PMBB" = "")) +
  coord_flip() +
  facet_wrap(~Hugo_symbol, nrow = 6, scale = "fixed") +
  expand_limits(y=0) +
  theme_Publication() +
  theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"))



# Create % patients with CPG PLP plot separated by source call, and include fractions as text 

perc_plot <- all_source_cts %>%
  filter(Hugo_symbol %in% sig_cpgs_both) %>%
  ggplot(aes(x = perc_plp, y = factor(cohort), fill = final_call_source, label = fraction)) +
  geom_bar(stat = "identity", color = "black",
           show.legend = TRUE) + 
  geom_text(x = 2.2, hjust = 0, size = 4, fontface = 2) +
  labs(x = "% Cohort P/LP", y = NULL, fill = NULL) + 
  scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                              "PMBB" = NULL)) +
  guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(~Hugo_symbol, nrow = 6, scale = "fixed") +
  expand_limits(x=2) +
  coord_cartesian(clip = 'off') +
  theme_Publication() +
  theme(plot.margin = unit(c(2,4,1,1), "lines"),
        legend.position = c(0.5, 1.07)) +
 # scale_fill_npg() +
  scale_fill_manual( values = c('#E64B35FF', '#4DBBD5FF'),
    limits = c('ClinVar', 'InterVar')
  )


# Merge plots and write to output

tiff(file.path(plot_dir, "sig-CPG-enrichment-PBTA-vs-control.tiff"),
     width = 10, height = 10, units = "in", res = 300)

ggarrange(pval_plot, enr_plot, perc_plot,
          nrow = 1, widths = c(2,1.25,1.5))

dev.off()



# Read in enrichment by CPG and histology files

hist_cpg_enr_gnomad <- read_tsv(cpg_hist_enr_gnomad_file)
hist_cpg_enr_pmbb <- read_tsv(cpg_hist_enr_pmbb_file)

# Create a 'dummy' data frame of enrichment results for CBTN to use as reference

hist_cpg_enr_pbta <- hist_cpg_enr_gnomad %>%
  select(plot_group, Hugo_symbol, n, n_without, hist_n, oddsRatio, pvalue, CI_lower, CI_upper, adjusted_p) %>%
  mutate(cohort = "PBTA",
         oddsRatio = NA,
         pvalue = NA,
         CI_lower = NA, 
         CI_upper = NA,
         adjusted_p = NA) %>%
  rename("n_plp" = n, "n_no_plp" = n_without)

# Subset gnomAD results for merging

hist_cpg_enr_gnomad <- hist_cpg_enr_gnomad %>%
  select(plot_group, Hugo_symbol, 
         count_with_plp_gnomAD, count_without_plp_gnomAD, 
         hist_n, oddsRatio, pvalue, CI_lower, CI_upper, adjusted_p) %>%
  mutate(cohort = "gnomAD") %>%
  rename("n_plp" = count_with_plp_gnomAD, 
         "n_no_plp" = count_without_plp_gnomAD)

# Subset PMBB results

hist_cpg_enr_pmbb <- hist_cpg_enr_pmbb %>%
  select(plot_group, Hugo_symbol, count_with_plp_PMBB_noTumor, 
         count_without_plp_PMBB_noTumor, 
         hist_n, oddsRatio, pvalue, CI_lower, CI_upper, adjusted_p) %>%
  mutate(cohort = "PMBB") %>%
  rename("n_plp" = count_with_plp_PMBB_noTumor, 
         "n_no_plp" = count_without_plp_PMBB_noTumor)

# merge all results and calculate % patients by histology with CPG PLP

hist_cpg_enr_all <- hist_cpg_enr_pbta %>%
  bind_rows(hist_cpg_enr_gnomad, hist_cpg_enr_pmbb) %>%
  mutate(plot_group = case_when(
    plot_group == "Atypical Teratoid Rhabdoid Tumor" ~ "ATRT",
    TRUE ~ plot_group
  )) %>%
  filter(!is.na(n_plp)) %>%
  mutate(perc_plp = n_plp/(n_plp + n_no_plp) * 100,
         hist_cpg = paste(plot_group, Hugo_symbol, sep = ":"),
         fraction = paste(round(n_plp,0), round(n_plp+n_no_plp,0), sep = "/"),
         cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA"))) %>%
  arrange(desc(-log10(pvalue))) %>%
  mutate(hist_cpg = factor(hist_cpg, unique(hist_cpg)))

# identify significantly enriched CPGs by histology relative to gnomAD and PMBB cohort, and retain those common to both lists

sig_hist_cpgs_gnomad <- hist_cpg_enr_all %>%
  filter(cohort == "gnomAD" & adjusted_p < 0.05 & adjusted_p > 0) %>%
  pull(hist_cpg)

sig_hist_cpgs_pmbb <- hist_cpg_enr_all %>%
  filter(cohort == "PMBB" & adjusted_p < 0.05 & adjusted_p > 0) %>%
  pull(hist_cpg)

sig_hist_cpgs_both <- intersect(sig_hist_cpgs_gnomad, sig_hist_cpgs_pmbb)



# Create data frame of PBTA CPG PLP counts by gene, histology, and call source 

pbta_source_hist_cts <- plp %>%
  rename("Kids_First_Biospecimen_ID_normal" = Kids_First_Biospecimen_ID) %>%
  left_join(hist[, c("Kids_First_Biospecimen_ID_normal", "plot_group")], by = "Kids_First_Biospecimen_ID_normal") %>%
  
  #Abbreviate ATRT for plots
  
  mutate(plot_group = case_when(
    plot_group == "Atypical Teratoid Rhabdoid Tumor" ~ "ATRT",
    TRUE ~ plot_group
  )) %>%
  
  # Classficy call source as ClinVar or Intervar
  
  mutate(final_call_source = case_when(grepl("ClinVar", Reasoning_for_call) ~ "ClinVar",
                                       grepl("InterVar", Reasoning_for_call) ~ "InterVar")) %>%
  
  # Get PLP counts by gene, plot group, and call source
  
  count(Hugo_Symbol, plot_group, final_call_source) %>%
  
  # add plot group n
  
  left_join(hist_cpg_enr_all[,c("plot_group", "hist_n")], by = c("plot_group")) %>%
  distinct(Hugo_Symbol, plot_group, final_call_source, .keep_all = T) %>%
  rename("n_plp" = n,
         "Hugo_symbol" = Hugo_Symbol) %>%
  mutate(n_no_plp = 838-hist_n,
         perc_plp = n_plp/hist_n*100,
         cohort = "PBTA")

# Merge with data from control cohorts

all_source_hist_cts <- hist_cpg_enr_all %>%
  filter(cohort != "PBTA") %>%
  
  # no call source or hist counts for controls 
  mutate(final_call_source = NA, hist_n = NA) %>%
  dplyr::select(Hugo_symbol, plot_group, cohort, final_call_source, n_plp, hist_n, n_no_plp, perc_plp) %>%
  bind_rows(pbta_source_hist_cts) %>%
  mutate(hist_cpg = paste(plot_group, Hugo_symbol, sep = ":")) %>%
  left_join(hist_cpg_enr_all[,c("cohort", "hist_cpg", "fraction")], by = c("cohort", "hist_cpg")) %>%
  distinct(hist_cpg, cohort, final_call_source, .keep_all = T) %>%
  mutate(cohort = factor(cohort, c("gnomAD", "PMBB", "PBTA")),
         hist_cpg = factor(hist_cpg, unique(hist_cpg_enr_all$hist_cpg)))
  


# create oddsRatio pvalue plot

hist_pval_plot <- hist_cpg_enr_all %>% 
  filter(hist_cpg %in% sig_hist_cpgs_both) %>%
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = pvalue, y = factor(cohort))) +
  geom_point(size = 3, show.legend = FALSE, 
             color = "#00A087FF") + 
  labs(x = "-log10(p)", y = "") + 
  geom_vline(xintercept = -log10(0.05/nrow(hist_cpg_enr_gnomad)), linetype = "dashed") + 
#  scale_y_discrete(labels=c("PBTA" = "", "gnomAD" = "",
#                            "PMBB" = "")) +
  facet_wrap(~hist_cpg, nrow = length(sig_hist_cpgs_both), scale = "fixed") +
  coord_cartesian(xlim = c(0,27)) + 
  theme_Publication() +
  theme(plot.margin = unit(c(3,1,1,0.5), "lines")) 


# create oddsRatio plot

hist_enr_plot <- hist_cpg_enr_all %>% 
  filter(hist_cpg %in% sig_hist_cpgs_both) %>%
  ggplot(aes(x = factor(cohort), y = oddsRatio)) +
  geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, 
                show.legend = FALSE, color = "#00A087FF") +
  labs(y = "Odds Ratio (95% CI)", x = NULL) + 
  scale_x_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                            "PMBB" = "")) +
  coord_flip() +
  facet_wrap(~hist_cpg, nrow = length(sig_hist_cpgs_both), scale = "fixed") +
  expand_limits(y=0) +
  theme_Publication() +
  theme(plot.margin = unit(c(3,0.5,1,1), "lines")) 


# Plot PLP frequency by gene and histology

hist_perc_plot <- all_source_hist_cts %>%
  filter(hist_cpg %in% sig_hist_cpgs_both) %>%
  ggplot(aes(x = perc_plp, y = factor(cohort), fill = final_call_source, label = fraction)) +
  geom_bar(stat = "identity", color = "black",
           show.legend = TRUE) + 
  
  geom_text(x = 70, hjust = 0, size = 4,
            fontface = 2) +
  labs(x = "% Cohort P/LP", y = NULL, fill = NULL) + 
  scale_y_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                            "PMBB" = "")) +
  guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(~hist_cpg, nrow = length(sig_hist_cpgs_both), scale = "fixed") +
  expand_limits(x=65) +
  coord_cartesian(clip = 'off') +
  theme_Publication() +
  theme(plot.margin = unit(c(3,4.5,1,1), "lines"),
        legend.position = c(0.5,1.035)) +
  scale_fill_manual( values = c('#E64B35FF', '#4DBBD5FF'),
                     limits = c('ClinVar', 'InterVar'))

# merge plots and save to output

png(file.path(plot_dir, "sig-hist-CPG-enrichment-PBTA-vs-control.png"),
     width = 11, height = 18, units = "in", res = 300)

ggarrange(hist_pval_plot, hist_enr_plot, hist_perc_plot,
          nrow = 1, widths = c(1.75,1.25,1.5))

dev.off()

