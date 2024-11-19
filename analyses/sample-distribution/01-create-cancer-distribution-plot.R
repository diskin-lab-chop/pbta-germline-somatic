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
plp_all <- read_tsv(file.path(root_dir, "analyses",
                             "bed-intersect", "results",
                             "pbta-merged-plp-variants-autogvp-abridged-wxs-exome-filtered-20bp_padded.tsv")) %>%
  filter(Kids_First_Biospecimen_ID_normal %in% hist$Kids_First_Biospecimen_ID_normal) %>%
  # determine whether final call was due to clinvar or intervar
  mutate(final_call_source = case_when(autogvp_call == "Likely_pathogenic" & autogvp_call_reason == "ClinVar" ~ "ClinVar - Likely Pathogenic",
                                       autogvp_call == "Pathogenic" & autogvp_call_reason == "ClinVar" ~ "ClinVar - Pathogenic",
                                       autogvp_call == "Likely_pathogenic" & autogvp_call_reason == "InterVar" ~ "InterVar - Likely Pathogenic",
                                       autogvp_call == "Pathogenic" & autogvp_call_reason == "InterVar" ~ "InterVar - Pathogenic"))

plp_cpg <- plp_all %>%
  filter(gene_symbol_vep %in% cpg)

# add n per group in label  
hist_counts <- hist %>%
  count(plot_group) %>%
  inner_join(hist) %>%
  mutate(plot_group_n = paste0(plot_group, " (n = ", n, ")")) %>%
  mutate(plot_group_n = fct_rev(fct_infreq(plot_group_n)),) %>%
  unique()

plot_group_palette <- hist_counts$plot_group_hex
names(plot_group_palette) <- hist_counts$plot_group_n

if(!interactive()) pdf(NULL)

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
  dplyr::select(final_call_source) %>%
  table(dnn = c("final_call_source", "n")) %>%
  as.data.frame() %>%
  dplyr::rename("final_call_source" = final_call_source.1) %>%
  mutate(final_call_source = fct_reorder(final_call_source, Freq))
plp_all_genes

plp_cpgs <- plp_cpg %>%
  dplyr::select(final_call_source) %>%
  table(dnn = c("final_call_source", "n")) %>%
  as.data.frame() %>%
  dplyr::rename("final_call_source" = final_call_source.1) %>%
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
  left_join(hist[,c("Kids_First_Biospecimen_ID_normal", "plot_group", "plot_group_hex")], by = "Kids_First_Biospecimen_ID_normal")

plp_sv <- read_tsv(file.path(data_dir, "pbta_germline_svs.tsv")) %>%
  left_join(hist[,c("Kids_First_Biospecimen_ID_normal", "plot_group", "plot_group_hex")], by = "Kids_First_Biospecimen_ID_normal") %>%
  dplyr::rename(gene_symbol_vep = Hugo_Symbol_cpg)

hist_snv_cpg <- plp_cpg %>%
  distinct(Kids_First_Biospecimen_ID_normal, .keep_all = T) %>%
  dplyr::filter(!Kids_First_Biospecimen_ID_normal %in% plp_sv$Kids_First_Biospecimen_ID_normal) %>%
  count(plot_group) 

hist_plp_cpg <- plp_sv %>%
  distinct(Kids_First_Biospecimen_ID_normal, .keep_all = T) %>%
  dplyr::mutate(plot_group = as.factor(plot_group)) %>%
  count(plot_group) %>%
  bind_rows(hist_snv_cpg) %>%
  group_by(plot_group) %>%
  summarize(n = sum(n)) %>%
  dplyr::rename('plp_cpg_n' = n) %>%
  full_join(hist_counts %>% distinct(plot_group, .keep_all = TRUE) %>%
              dplyr::select(plot_group, n,
                            plot_group_n)) %>%
  arrange(plot_group) %>%
  dplyr::mutate(plp_cpg_n = case_when(
    is.na(plp_cpg_n) ~ 0,
    TRUE ~ plp_cpg_n
  )) %>%
  distinct(plot_group, .keep_all = TRUE) %>%
  mutate(freq = plp_cpg_n/n*100) %>%
  dplyr::mutate(plot_group_n = glue::glue("{plot_group} (n = {plp_cpg_n}/{n})"))

# rename plot group palette colors to match `plot_group_n` in `hist_plp_cpg`
plot_group_palette <- unique(plot_group_palette)
names(plot_group_palette) <- hist_plp_cpg$plot_group_n

# create count and freq plots by group

tiff(file.path(plot_dir, "histology-CPG-PLP-count.tiff"), height = 1800, width = 2800, res = 300)
hist_plp_cpg %>%
  arrange(plp_cpg_n) %>%
  ggplot(aes(x = factor(plot_group_n, plot_group_n), y = plp_cpg_n, fill = factor(plot_group_n))) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) + 
  labs(x = "", y = "Number of patients with germline CPG P/LP variant") +
  scale_fill_manual(values = plot_group_palette) +
  coord_flip() + 
  theme_Publication()
dev.off()

tiff(file.path(plot_dir, "histology-CPG-PLP-freq.tiff"), height = 1800, width = 2500, res = 300)
hist_plp_cpg %>%
  arrange(freq) %>%
  ggplot(aes(x = factor(plot_group_n, plot_group_n), y = freq, fill = factor(plot_group_n))) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) + 
  labs(x = "", y = "% Patients with germline CPG P/LP variant") +
  scale_fill_manual(values = plot_group_palette) +
  coord_flip() + 
  theme_Publication()
dev.off()

# obtain gene-level count and freq of PLP variants by plot group 

hist_gene_plp_cpg <- plp_cpg %>%
  distinct(Kids_First_Biospecimen_ID_normal, gene_symbol_vep, .keep_all = T) %>%
  count(plot_group, gene_symbol_vep) %>%
  dplyr::rename('plp_cpg_n' = n) %>%
  filter(!is.na(plot_group)) %>%
  left_join(unique(hist_counts[,c("plot_group", "n", "plot_group_n")], by = "plot_group")) %>%
  mutate(freq = plp_cpg_n/n) %>%
  arrange(plot_group, freq) %>%
  mutate(gene_symbol_vep = factor(gene_symbol_vep, unique(gene_symbol_vep)))

hist_gene_snv_cpg <- plp_cpg %>%
  count(plot_group, gene_symbol_vep) 

hist_gene_plp_cpg <- plp_sv %>%
  count(plot_group, gene_symbol_vep) %>%
  bind_rows(hist_gene_snv_cpg) %>%
  group_by(plot_group, gene_symbol_vep) %>%
  summarize(n = sum(n)) %>%
  dplyr::rename('plp_cpg_n' = n) %>%
  left_join(unique(hist_counts[,c("plot_group", "n", "plot_group_n")], by = "plot_group")) %>%
  mutate(freq = plp_cpg_n/n*100)

# Create gene-by-histology freq plot

png(file.path(plot_dir, "CPG-PLP-freq-by-histology.png"), height = 5800, width = 5800, res = 300)
hist_gene_plp_cpg %>%
  mutate(plot_group_n = factor(plot_group_n, unique(hist_gene_plp_cpg$plot_group_n)),
         gene_symbol_vep = reorder_within(gene_symbol_vep, freq, plot_group_n),
         plot_group_n = case_when(
           plot_group_n == "Atypical Teratoid Rhabdoid Tumor (n = 27)" ~ "ATRT (n = 27)",
           TRUE ~ plot_group_n
         )) %>%
  ggplot(aes(x = gene_symbol_vep, y = freq)) +
  geom_point(size = 3, show.legend = FALSE) + 
  geom_segment(aes(x=gene_symbol_vep, xend=gene_symbol_vep, y=0, yend=freq),
               linewidth = 1,
               show.legend = FALSE) +
  labs(x = "", y = "Proportion of patients with germline PLP variant") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~plot_group_n, scale = "free") +
  theme_Publication()
dev.off()


