library(ggplot2)
library(tidyverse) 
library(tidytext)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "sample-distribution")
plot_dir <- file.path(analysis_dir, "plots/")
input_dir <- file.path(analysis_dir, "input")

# If the results directory does not exist, create it
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# source publication theme
source(file.path(root_dir, "figures", "theme.R"))

# read in histologies file
hist <- read_tsv(file.path(root_dir, 
                           "analyses", 
                           "collapse-tumor-histologies", 
                           "results",
                           "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv"))

# read in CPG list
cpg <- read_lines(file.path(root_dir, "analyses", "oncokb-annotation", "input", "cpg.txt"))


# read in plp file
plp_all <- read_tsv(file.path(data_dir, "pbta_germline_plp_calls.tsv")) %>%
  filter(Kids_First_Biospecimen_ID %in% hist$Kids_First_Biospecimen_ID_normal) %>%
  # determine whether final call was due to clinvar or intervar
  mutate(final_call_source = case_when(grepl("ClinVar Likely pathogenic", Reasoning_for_call) ~ "ClinVar - Likely Pathogenic",
                                       grepl("ClinVar Pathogenic", Reasoning_for_call) ~ "ClinVar - Pathogenic",
                                       grepl("InterVar_Recalculated Likely_pathogenic", Reasoning_for_call) ~ "InterVar - Likely Pathogenic",
                                       grepl("InterVar_Recalculated Pathogenic", Reasoning_for_call) ~ "InterVar - Pathogenic"))

plp_cpg <- plp_all %>%
  filter(Hugo_Symbol %in% cpg)


# add n per group in label  
hist_counts <- hist %>%
  count(plot_group) %>%
  inner_join(hist) %>%
  mutate(plot_group_n = paste0(plot_group, " (n = ", n, ")")) %>%
  mutate(plot_group_n = fct_rev(fct_infreq(plot_group_n)),) %>%
  unique()

plot_group_palette <- hist_counts$plot_group_hex
names(plot_group_palette) <- hist_counts$plot_group_n

tiff(file.path(plot_dir, "histology-distribution.tiff"), height = 1500, width = 2500, res = 300)
ggplot(hist_counts, aes(x = plot_group_n, fill = plot_group_n)) +
  geom_bar(color = "black", show.legend = FALSE) +
  scale_fill_manual(values = plot_group_palette) + 
  xlab("Histology") +
  ylab("Number of patients with tumor diagnosis") +
  coord_flip() + 
  theme_Publication()
dev.off()


# current + x01 cohort 
x01 <- read_tsv(file.path(input_dir, "838_OpenPBTAcohort_1720_X01cohort.tsv")) %>%
  filter(!is.na(plot_group))


# add n per group in label  
x01_counts <- x01 %>%
  count(plot_group, cohort) %>%
  group_by(plot_group) %>%
  # make sure OpenPBTA is first - relevel factors after this for the plot
  summarise(n = str_c(n, sep = ", ", collapse = ", ")) %>%
  inner_join(x01) %>%
  mutate(plot_group_n = paste0(plot_group, " (n = ", n, ")"),
         # reorder plot groups in descending order
         plot_group_n = fct_rev(fct_infreq(plot_group_n)),
        # now relevel cohort to make OpenPBTA appear first
        cohort = fct_relevel(cohort, "PBTA X01", "OpenPBTA"))


tiff(file.path(plot_dir, "x01-histology-distribution.tiff"), height = 1800, width = 3600, res = 300)
ggplot(x01_counts, aes(x = plot_group_n, fill = cohort)) +
  geom_bar(color = "black", show.legend = TRUE) +
  scale_fill_manual(values = c("OpenPBTA" = "black", "PBTA X01" = "grey")) + 
  xlab("Histology") +
  ylab("Number of patients with tumor diagnosis") +
  coord_flip() + 
  theme_Publication()
dev.off()

# summary of variants for all genes or just CPGs
plp_all_genes <- plp_all %>%
  select(final_call_source) %>%
  table(dnn = c("final_call_source", "n")) %>%
  as.data.frame() %>%
  rename("final_call_source" = final_call_source.1) %>%
  mutate(final_call_source = fct_reorder(final_call_source, Freq))
plp_all_genes

plp_cpgs <- plp_cpg %>%
  select(final_call_source) %>%
  table(dnn = c("final_call_source", "n")) %>%
  as.data.frame() %>%
  rename("final_call_source" = final_call_source.1) %>%
  mutate(final_call_source = fct_reorder(final_call_source, Freq))
plp_cpgs

# create plots
dfs <- list("plp_all_genes" = plp_all_genes, "plp_cpgs" = plp_cpgs)

for (df in names(dfs)) {
  
  # add titles to print
  if (df == "plp_all_genes"){
    title <- "All genes"
  }
  
  if (df == "plp_cpgs"){
    title <- "CPGs"
  }
    
  dev.set(dev.next())
  tiff(file.path(paste0(plot_dir, df, "-calls.tiff")), height = 1000, width = 1600, res = 300)
  
  print(
    dfs[[df]] %>%
    arrange(Freq) %>%
    mutate(group=factor(final_call_source, final_call_source)) %>%
    ggplot( aes(x=final_call_source, y=Freq) ) +
    geom_segment( aes(x=final_call_source ,xend=final_call_source, y=0, yend=Freq), color="grey") +
    geom_point(size=3, color="black") +
    coord_flip() +
    theme_Publication() +
    theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
    ) +
    ylab("Number of germline variants") +
    xlab("Final call") +
    ggtitle(title) +
    #left align
    theme(plot.title = element_text(hjust = 0))   
  )
  dev.off()
  
}


# append plot group to plp_cpg
plp_cpg <- plp_cpg %>%
  rename("Kids_First_Biospecimen_ID_normal" = Kids_First_Biospecimen_ID) %>%
  left_join(hist[,c("Kids_First_Biospecimen_ID_normal", "plot_group", "plot_group_hex")], by = "Kids_First_Biospecimen_ID_normal")

# obtain counts and freq of CPG PLP variants by plot group 

hist_plp_cpg <- plp_cpg %>%
  distinct(Kids_First_Biospecimen_ID_normal, .keep_all = TRUE) %>%
  count(plot_group) %>%
  dplyr::rename('plp_cpg_n' = n) %>%
  filter(!is.na(plot_group)) %>%
  left_join(hist_counts[,c("plot_group", "n", "plot_group_n")], by = "plot_group") %>%
  distinct(plot_group, .keep_all = TRUE) %>%
  mutate(freq = plp_cpg_n/n*100)

# create count and freq plots by group

tiff(file.path(plot_dir, "histology-CPG-PLP-count.tiff"), height = 1800, width = 2800, res = 300)
hist_plp_cpg %>%
  arrange(plp_cpg_n) %>%
  ggplot(aes(x = factor(plot_group_n, plot_group_n), y = plp_cpg_n, fill = factor(plot_group_n))) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) + 
  labs(x = "", y = "Number of patients with germline CPG P-LP variant") +
  scale_fill_manual(values = plot_group_palette) +
  coord_flip() + 
  theme_Publication()
dev.off()

tiff(file.path(plot_dir, "histology-CPG-PLP-freq.tiff"), height = 1800, width = 2500, res = 300)
hist_plp_cpg %>%
  arrange(freq) %>%
  ggplot(aes(x = factor(plot_group_n, plot_group_n), y = freq, fill = factor(plot_group_n))) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) + 
  labs(x = "", y = "% Patients with germline CPG P-LP variant") +
  scale_fill_manual(values = plot_group_palette) +
  coord_flip() + 
  theme_Publication()
dev.off()

# obtain gene-level count and freq of PLP variants by plot group 

hist_gene_plp_cpg <- plp_cpg %>%
  count(plot_group, Hugo_Symbol) %>%
  dplyr::rename('plp_cpg_n' = n) %>%
  filter(!is.na(plot_group)) %>%
  left_join(hist_counts[,c("plot_group", "n", "plot_group_n")], by = "plot_group") %>%
  mutate(freq = plp_cpg_n/n) %>%
  arrange(plot_group, freq) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, unique(Hugo_Symbol)))

# Create freq plot

png(file.path(plot_dir, "CPG-PLP-freq-by-histology.png"), height = 5700, width = 5500, res = 300)
hist_gene_plp_cpg %>%
  mutate(plot_group_n = factor(plot_group_n, unique(hist_gene_plp_cpg$plot_group_n)),
         Hugo_Symbol = reorder_within(Hugo_Symbol, freq, plot_group_n)) %>%
  ggplot(aes(x = Hugo_Symbol, y = freq)) +
  geom_point(size = 3, show.legend = FALSE) + 
  geom_segment(aes(x=Hugo_Symbol, xend=Hugo_Symbol, y=0, yend=freq),
               linewidth = 1,
               show.legend = FALSE) +
  labs(x = "", y = "Proportion of patients with germline PLP variant") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~plot_group_n, scale = "free") +
  theme_Publication()
dev.off()

