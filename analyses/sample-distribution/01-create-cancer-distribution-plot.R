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

# source publication theme
source(file.path(root_dir, "figures", "theme.R"))
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
as.data.frame(table(hist_counts$plot_group_n))

plot_group_palette <- hist_counts$plot_group_hex
names(plot_group_palette) <- hist_counts$plot_group_n

"#8f8fbf"

tiff(file.path(plot_dir, "histology-distribution.tiff"), height = 2000, width = 3000, res = 300)
ggplot(hist_counts, aes(x = plot_group_n, fill = plot_group_n)) +
  geom_bar(color = "black", show.legend = FALSE) +
  scale_fill_manual(values = plot_group_palette) + 
  xlab("Histology") +
  ylab("Number of patients with tumor diagnosis") +
  coord_flip() + 
  theme_Publication()
dev.off()

