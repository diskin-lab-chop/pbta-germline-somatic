library(ggplot2)
library(tidyverse) 

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "variant-distribution")
plot_dir <- file.path(analysis_dir, "plots")

# If the results directory does not exist, create it
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# source publication theme
source(file.path(root_dir, "figures", "theme.R"))

# read in cpg list
cpg <- read_lines(file.path(root_dir, "analyses", "oncokb-annotation", "input", "cpg.txt"))

# read in histology file
hist <- read_tsv(file.path(root_dir, 
                           "analyses", 
                           "collapse-tumor-histologies", 
                           "results",
                           "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv"))

# merge with germline plp file
cpg_plp <- read_tsv(file.path(data_dir, "pbta_germline_plp_calls.tsv")) %>%
  filter(Kids_First_Biospecimen_ID %in% hist$Kids_First_Biospecimen_ID_normal) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID) %>%
  # unique one gene per patient if multiple variants (in this case we have 1 normal per pt, so can use bs_id)
  select(Kids_First_Biospecimen_ID_normal, Hugo_Symbol) %>%
  unique() %>%
  # filter for CPGs
  filter(Hugo_Symbol %in% cpg)

# combine plp var with hist data 
cpg_plp_plot_gp <- cpg_plp %>%
  left_join(hist) 

# get N per histology per gene
n_per_hist <- cpg_plp_plot_gp %>%
  count(Hugo_Symbol, plot_group) %>%
  dplyr::rename(hist_n = n)

# create N per gene desc order for plot
n_per_gene <- cpg_plp %>%
  left_join(hist) %>%
  count(Hugo_Symbol) %>%
  dplyr::rename(gene_n = n) %>%
  arrange(desc(gene_n)) %>%
  # create plot order desc by gene
  mutate(plot_order = row_number()) %>%
  # add back histology info
  left_join(n_per_hist) %>%
  left_join(cpg_plp_plot_gp[,c("plot_group", "plot_group_hex")]) %>%
  distinct()

# palette
plot_group_palette <- n_per_gene$plot_group_hex
names(plot_group_palette) <- n_per_gene$plot_group


tiff(file.path(plot_dir, "variant-distribution-by-histology.tiff"), height = 1800, width = 4000, res = 300)
ggplot(data = n_per_gene, aes(x = reorder(Hugo_Symbol, plot_order), y = hist_n, fill = plot_group)) +
  geom_col(col = "black", size = 0.4) +
  scale_fill_manual(values = plot_group_palette) + 
  ylab("Number of patients with P/LP variants") +
  xlab("Gene") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0,20) +
  labs(fill = "Tumor diagnosis")
dev.off()








