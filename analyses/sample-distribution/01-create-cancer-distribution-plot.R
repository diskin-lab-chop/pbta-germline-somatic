library(ggplot2)
library(tidyverse) 

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-distribution")
plot_dir <- file.path(analysis_dir, "plots")

# If the results directory does not exist, create it
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}


#source("~/Box Sync/D3B-share/manuscripts/02_accepted/Lilly_CBTN_infrastructure/figures/theme.R")
hist <- read_tsv(file.path(root_dir, 
                           "analyses", 
                           "collapse-tumor-histologies", 
                           "results",
                           "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv"))
  

hist_counts <- hist %>%
  count(plot_group) %>%
  inner_join(hist) %>%
  mutate(plot_group_n = paste0(plot_group, " (n = ", n, ")")) %>%
  mutate(plot_group_n = fct_rev(fct_infreq(plot_group_n))) %>%
  unique()


plot_group_palette <- hist_counts$plot_group_hex
names(plot_group_palette) <- hist_counts$plot_group_n


tiff("~/Box Sync/D3B-share/manuscripts/02_accepted/Lilly_CBTN_infrastructure/figures/tumors_by_histology_revision.tiff", height = 2000, width = 3000, res = 300)
ggplot(tumor_counts, aes(x = cancer_group_display_n, fill = cancer_group_display_n)) +
  geom_bar(color = "black", show.legend = FALSE) +
  scale_fill_manual(values = cancer_group_palette) + 
  xlab("Histology") +
  ylab("Tumor count") +
  ylim(c(0,150)) +
  coord_flip() + 
  ggtitle("Number of tumors profiled in OpenPBTA") +
  theme_Publication()
dev.off()

