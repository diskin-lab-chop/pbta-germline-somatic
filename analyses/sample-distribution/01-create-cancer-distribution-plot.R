library(ggplot2)
library(tidyverse) 

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
  table(dnn = c("call", "n")) %>%
  as.data.frame() %>%
  mutate(final_call_source = fct_reorder(final_call_source, Freq))
plp_all_genes

plp_cpgs <- plp_cpg %>%
  select(final_call_source) %>%
  table(dnn = c("call", "n")) %>%
  as.data.frame() %>%
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


