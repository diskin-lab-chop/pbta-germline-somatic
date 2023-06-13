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
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

cpg_enr_file <- file.path(input_dir, "PBTA-all-CPG-PLP-enrichment.tsv")

cpg_enr <- read_tsv(cpg_enr_file) %>%
  dplyr::filter(cohort != "PMBB_noCancer") %>%
  dplyr::mutate(cohort = fct_relevel(cohort, 
                                     c("gnomAD", "PMBB_noTumor", "PBTA")))



# Create CPG enrichment FDR plot  

pval_plot <- cpg_enr %>% 
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = pvalue, y = factor(cohort))) +
  geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
  labs(x = "-log10(p)", y = "") + 
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") + 
  #coord_cartesian(xlim = c(0,24)) + 
  theme_Publication() +
  theme(plot.margin = unit(c(2,1,1,0), "lines"))


# Create CPG Odds Ratio plot 

enr_plot <- cpg_enr %>% 
  ggplot(aes(x = factor(cohort), y = OddsRatio)) +
  geom_point(size = 3, color = "#00A087FF",
             show.legend = FALSE) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, 
                show.legend = FALSE, color = "#00A087FF") +
  labs(y = "log10-Odds Ratio (95% CI)", x = NULL) + 
  scale_x_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                            "PMBB_noTumor" = "")) +
  coord_flip() +
  expand_limits(y=0) +
  theme_Publication() +
  theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"))



# Create % patients with CPG PLP plot separated by source call, and include fractions as text 

perc_plot <- cpg_enr %>%
  ggplot(aes(x = perc_plp, y = factor(cohort), label = fraction)) +
  geom_bar(stat = "identity", color = "black",
           show.legend = TRUE) + 
  geom_text(x = 22, hjust = 0, size = 4, fontface = 2) +
  labs(x = "% Cohort P/LP", y = NULL, fill = NULL) + 
  scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                            "PMBB_noTumor" = NULL)) +
  guides(fill = guide_legend(nrow = 1)) +
  expand_limits(x=2) +
  coord_cartesian(clip = 'off') +
  theme_Publication() +
  theme(plot.margin = unit(c(2,5,1,1), "lines"))


# Merge plots and write to output

pdf(file.path(plot_dir, "sig-CPG-enrichment-PBTA-vs-control.tiff"),
     width = 8, height = 2.5)

ggarrange(pval_plot, enr_plot, perc_plot,
          nrow = 1, widths = c(2,1.25,1.5))

dev.off()
