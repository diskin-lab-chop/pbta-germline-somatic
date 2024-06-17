################################################################################
# 02-cohort-summary-plots.R
# Generate cohort summary plots for PBTA germline project
# Author: Ryan Corbett
# Usage: Rscript 02-cohort-summary-plots.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(colorblindr)
  library(ggpubr)
})

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "collapse-tumor-histologies")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
plot_dir <- file.path(analysis_dir, "plots")

# If the results directory does not exist, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

source(file.path(root_dir, "figures", "theme.R"))

# set relevant file paths
hist_file <- file.path(results_dir, "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

opc_hist_file <- file.path(data_dir, "histologies.tsv")

rna_specimens_file <- file.path(data_dir, 
                                "independent-specimens.rnaseqpanel.primary-plus.tsv")

methylation_specimens_file <- file.path(data_dir,
                                        "independent-specimens.methyl.primary-plus.tsv")

cptac_proteomics_file <- file.path(data_dir,
                                  "cptac-protein-imputed-prot-expression-abundance.tsv.gz")

hope_proteomics_file <- file.path(data_dir,
                                  "hope-protein-imputed-prot-expression-abundance.tsv.gz")

# Load OPC histologies file
opc_hist <- read_tsv(opc_hist_file)

# Define RNA, methylation, and proteomics specimens
rna_specimens <- opc_hist %>%
  dplyr::filter(experimental_strategy == "RNA-Seq") %>%
  distinct(match_id, .keep_all = TRUE) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID)

methylation_specimens <- opc_hist %>%
  dplyr::filter(experimental_strategy == "Methylation") %>%
  distinct(match_id, .keep_all = TRUE) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_methyl = Kids_First_Biospecimen_ID)

proteo_specimens <- opc_hist %>%
  dplyr::filter(experimental_strategy == "Whole Cell Proteomics") %>%
  dplyr::select(match_id, Kids_First_Biospecimen_ID) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_proteo = Kids_First_Biospecimen_ID) %>%
  distinct(match_id, .keep_all = TRUE)
  
# Load cohort hist file and append RNA, methylation, and proteomics IDs when available
hist <- read_tsv(hist_file) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_tumor) %>%
  left_join(opc_hist %>%
              dplyr::select(Kids_First_Biospecimen_ID, match_id)) %>%
  left_join(rna_specimens %>%
              dplyr::select(Kids_First_Biospecimen_ID_RNA, match_id)) %>%
  left_join(methylation_specimens %>%
              dplyr::select(Kids_First_Biospecimen_ID_methyl, match_id)) %>%
  left_join(proteo_specimens) %>%
  dplyr::mutate(plot_group = case_when(
    plot_group == "Atypical Teratoid Rhabdoid Tumor" ~ "AT/RT",
    TRUE ~ plot_group
  )) %>%
  dplyr::mutate(RNA = case_when(
    !is.na(Kids_First_Biospecimen_ID_RNA) ~ "Yes",
    TRUE ~ "No"
  )) %>%
  dplyr::mutate(Methylation = case_when(
    !is.na(Kids_First_Biospecimen_ID_methyl) ~ "Yes ",
    TRUE ~ "No"
  )) %>%
  dplyr::mutate(Proteomics = case_when(
    !is.na(Kids_First_Biospecimen_ID_proteo) ~ "Yes  ",
    TRUE ~ "No"
  ))
  
# Get participant count by plot group
plot_group_cts <- hist %>%
  count(plot_group) %>%
  dplyr::mutate(group_n = glue::glue("{plot_group} (n={n})")) %>%
  left_join(hist %>% dplyr::select(plot_group, plot_group_hex)) %>%
  distinct()

# create circos plot df with relevant data columns
plot_df <- hist %>%
  left_join(plot_group_cts %>% dplyr::select(-n, -plot_group_hex)) %>%
  dplyr::select(Kids_First_Participant_ID,
                group_n,
                RNA, Proteomics, Methylation) %>%
  column_to_rownames("Kids_First_Participant_ID") %>%
  dplyr::arrange(group_n,
                 RNA, Proteomics, Methylation)

# define color palettes
plot_group_cols <- plot_group_cts$plot_group_hex
names(plot_group_cols) <- plot_group_cts$group_n

rna_cols <- c("gray80", "whitesmoke")
names(rna_cols) <- c("Yes", "No")

prot_cols <- c("gray40", "whitesmoke")
names(prot_cols) <- c("Yes  ", "No")

methyl_cols <- c("gray10", "whitesmoke")
names(methyl_cols) <- c("Yes ", "No")

cols_all <- c(plot_group_cols,
              rna_cols, prot_cols, methyl_cols)

# split circos plot by plot group
split <- factor(plot_df$group_n)

pdf(NULL)

# Generate circos plot
pdf(file.path(plot_dir, "circos-plot.pdf"), width = 7, height = 7)

circos.clear()
circos.heatmap(plot_df,
               split=split,
               col = cols_all)

# add border colors to sectors
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = 1)
  circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
              CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
              col = NA, border = "black")
}

## add legends
lgd_plot_group = Legend(title = "Histology", 
                        at =  sort(names(plot_group_cols)), 
                        legend_gp = gpar(fill = plot_group_cols[sort(names(plot_group_cols))]))

lgd_rna = Legend(title = "RNA-seq", 
                    at =  names(rna_cols), 
                    legend_gp = gpar(fill = rna_cols))

lgd_prot = Legend(title = "Proteomics", 
                 at =  names(prot_cols), 
                 legend_gp = gpar(fill = prot_cols))

lgd_methyl = Legend(title = "DNA methylation", 
                 at =  names(methyl_cols), 
                 legend_gp = gpar(fill = methyl_cols))

## Placement of legend based on device size
h = dev.size()[2]
circle_size = unit(1, "snpc")

# Create legends
lgd_list1 = packLegend(lgd_plot_group, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_rna, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_prot, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list4 = packLegend(lgd_methyl, max_height = unit(0.9*h, "inch"), direction = "horizontal")

# Add legends on plot
draw(lgd_list1, x = unit(76, "mm"), y = unit(90, "mm"))
draw(lgd_list2, x = unit(118, "mm"), y = unit(112, "mm"))
draw(lgd_list3, x = unit(119, "mm"), y = unit(97, "mm"))
draw(lgd_list4, x = unit(124, "mm"), y = unit(82, "mm"))

circos.clear()
dev.off()

# Define df of sex estimate for pie chart
sex_df <- data.frame(labels = c(glue::glue("Female (n={sum(hist$germline_sex_estimate == 'Female')})"),
                                glue::glue("Male (n={sum(hist$germline_sex_estimate == 'Male')})"),
                                glue::glue("Unknown (n={sum(hist$germline_sex_estimate == 'Unknown')})")),
                 sizes = c(sum(hist$germline_sex_estimate == "Female"), 
                           sum(hist$germline_sex_estimate == "Male"),
                           sum(hist$germline_sex_estimate == "Unknown")))

# plot pie chart
sex_pie_chart <- ggplot(sex_df, aes(x = "", y = sizes, fill = labels)) +
  geom_bar(stat = "identity", width = 2,
           colour = "black",
           show.legend = TRUE) +
  coord_polar(theta = "y") +
  labs(fill = "Sex") +
  scale_fill_manual(values = c("#7CAE00", "#F8766D", "gray")) +
  theme_void() +
  theme(legend.text = element_text(size=22),
        legend.title = element_text(face = "bold",
                                    size = 22),
        legend.text.align = 0)

# Define ancestry unknown for pie chart
hist <- hist %>%
  dplyr::mutate(predicted_ancestry = case_when(
    is.na(predicted_ancestry) ~ "Unknown",
    TRUE ~ predicted_ancestry
  ))

# generate genetic ancestry pie chart df
anc_df <- data.frame(labels = c(glue::glue("AFR (n={sum(hist$predicted_ancestry == 'AFR')})"),
                                glue::glue("AMR (n={sum(hist$predicted_ancestry == 'AMR')})"),
                                glue::glue("EAS (n={sum(hist$predicted_ancestry == 'EAS')})"),
                                glue::glue("EUR (n={sum(hist$predicted_ancestry == 'EUR')})"),
                                glue::glue("SAS (n={sum(hist$predicted_ancestry == 'SAS')})"),
                                glue::glue("Unknown (n={sum(hist$predicted_ancestry == 'Unknown')})")),
                     
                     sizes = c(sum(hist$predicted_ancestry == "AFR"), 
                               sum(hist$predicted_ancestry == "AMR"),
                               sum(hist$predicted_ancestry == "EAS"),
                               sum(hist$predicted_ancestry == "EUR"),
                               sum(hist$predicted_ancestry == "SAS"),
                               sum(hist$predicted_ancestry == "Unknown")))

anc_pie_chart <- ggplot(anc_df, aes(x = "", y = sizes, fill = labels)) +
  geom_bar(stat = "identity", width = 2,
           colour = "black",
           show.legend = TRUE) +
  coord_polar(theta = "y") +
  labs(fill = "Predicted Ancestry") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73",
                               "#0072B2", "#D55E00", "gray")) +
  theme_void() +
  theme(legend.text = element_text(size=20),
        legend.title = element_text(face = "bold",
                                    size = 22),
        legend.box.just = "left")

# Generate tumor event pie chart
tumor_df <- data.frame(labels = c(glue::glue("Initial CNS tumor (n={sum(hist$tumor_descriptor == 'Initial CNS Tumor')})"),
                                glue::glue("Recurrence (n={sum(hist$tumor_descriptor == 'Recurrence')})"),
                                glue::glue("Deceased (n={sum(hist$tumor_descriptor == 'Deceased')})"),
                                glue::glue("Second Malignancy (n={sum(hist$tumor_descriptor == 'Second Malignancy')})"),
                                glue::glue("Progressive (n={sum(hist$tumor_descriptor == 'Progressive')})"),
                                glue::glue("Unknown (n={sum(hist$tumor_descriptor == 'Unknown')})")),
                     
                     sizes = c(sum(hist$tumor_descriptor == "Initial CNS Tumor"), 
                               sum(hist$tumor_descriptor == "Recurrence"),
                               sum(hist$tumor_descriptor == "Deceased"),
                               sum(hist$tumor_descriptor == "Second Malignancy"),
                               sum(hist$tumor_descriptor == "Progressive"),
                               sum(hist$tumor_descriptor == "Unknown")))

tumor_pie_chart <- ggplot(tumor_df, aes(x = "", y = sizes, fill = labels)) +
  geom_bar(stat = "identity", width = 2,
           colour = "black",
           show.legend = TRUE) +
  coord_polar(theta = "y") +
  labs(fill = "Tumor Event") +
  scale_fill_manual(values = c("#CC6677", "#DDCC77", "#117733",
                               "#332288", "#661100", "gray")) +
  theme_void() +
  theme(legend.text = element_text(size=16),
        legend.title = element_text(face = "bold",
                                    size = 22),
        legend.position = "right",
        legend.text.align = 0)

# Merge pie charts and save
all_pie <- ggarrange(sex_pie_chart,
                 anc_pie_chart,
                 tumor_pie_chart,
                 ncol = 1, align = "hv")

all_pie

ggsave(file.path(plot_dir, "demo-piecharts.pdf"),
       width = 6, height = 9, units = "in")

