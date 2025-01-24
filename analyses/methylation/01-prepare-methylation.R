library(tidyverse) 
library(biomaRt)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "methylation")
plot_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")


# set file paths

methyl_file <- file.path(data_dir,
                         "methyl-beta-values.rds")

methyl_annot_file <- file.path(data_dir,
                               "infinium.gencode.v39.probe.annotations.tsv.gz")

hist_file <- file.path(root_dir, "analyses",
                       "collapse-tumor-histologies", "results",
                       "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

opc_hist_file <- file.path(data_dir, "histologies.tsv")

cpg_file <- file.path(root_dir, "analyses",
                      "oncokb-annotation", "input",
                      "cpg.txt")

plp_file <- file.path(data_dir, 
                      "pbta-merged-plp-variants-autogvp-abridged.tsv")

plp_lowVAF_file <- file.path(data_dir,
                             "pbta-merged-plp-variants-autogvp-abridged-lowVAF-cpgs.tsv")

plp_sv_file <- file.path(data_dir, 
                         "pbta_germline_svs.tsv")


# First, check that methylation data files exist, and download if not

if (!file.exists(methyl_file) | !file.exists(methyl_annot_file)){
  
  system(glue::glue("bash {root_dir}/scripts/download-methyl.sh"))
  
}

# Wrangle histology data

opc_hist <- read_tsv(opc_hist_file, guess_max = 10000)

hist <- read_tsv(hist_file, guess_max = 10000) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_tumor) %>%
  # add match id
  left_join(opc_hist %>% dplyr::select(Kids_First_Biospecimen_ID,
                                       match_id),
            by = "Kids_First_Biospecimen_ID") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID) %>%
  # add matched methylation BS ID
  left_join(opc_hist %>% 
              dplyr::filter(experimental_strategy == "Methylation") %>%
              dplyr::select(match_id, Kids_First_Biospecimen_ID),
            by = "match_id") %>%
  distinct(Kids_First_Biospecimen_ID_tumor, .keep_all = TRUE) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_methyl = Kids_First_Biospecimen_ID) %>%
  write_tsv(file.path(results_dir,
                      "germline-primary-plus-tumor-histologies-methylation.tsv"))

# extract methylation BS IDs
methyl_ids <- hist %>%
  dplyr::filter(!is.na(Kids_First_Biospecimen_ID_methyl)) %>%
  pull(Kids_First_Biospecimen_ID_methyl)

# Load methylation data
methyl <- readRDS(methyl_file)

# filter for samples in cohort and probes with data
methyl <- methyl[, colnames(methyl) %in% c("Probe_ID", methyl_ids)]
methyl <- methyl[rowSums(!is.na(methyl[,-1])) != 0,]

# Obtain sample mean methylation beta values
num_iterations <- ceiling(ncol(methyl) / 100)

mean_beta <- rep(0, ncol(methyl))

# loop through columns in groups of 100 and calculate colMeans
for (i in 1:num_iterations) {
  # Determine the indices for the current iteration
  if (i == 1) {
    
    start_index <- (i - 1) * 100 + 2
    end_index <- min(i * 100, ncol(methyl))
    
  } else {
    
    start_index <- (i - 1) * 100 + 1
    end_index <- min(i * 100, ncol(methyl))
    
  }
  
  # Extract the elements for the current iteration
  mean_beta[(start_index):(end_index)] <- colMeans(methyl[,start_index:end_index], na.rm = T)

}

mean_beta_df <- data.frame(Kids_First_Biospecimen_ID_methyl = colnames(methyl)[-1],
                           mean_beta_value = mean_beta[-1])

write_tsv(mean_beta_df,
          file.path(results_dir, 
                    "pbta-germline-mean-sample-methyl-beta.tsv"))

# Load annotation file and filter for CPG-annotated probes 

cpgs <- read_lines(cpg_file)

annot <- read_tsv(methyl_annot_file) %>%
  # filter for probes in `methyl` and for probes in CPGs
  dplyr::filter(`Probe_ID` %in% methyl$Probe_ID,
                Gene_symbol %in% cpgs)

transcript_ids <- annot %>%
  dplyr::filter(!is.na(transcript_id)) %>%
  distinct(transcript_id) %>%
  pull(transcript_id)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get canonical yes/no for all transcript IDs
bm <- getBM(
  attributes = c("ensembl_transcript_id", "transcript_is_canonical"),
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl
)

bm <- bm %>%
  filter(transcript_is_canonical == 1)

# add canonical yes/no to annot df
annot <- annot %>%
  dplyr::mutate(canonical = case_when(
    transcript_id %in% bm$ensembl_transcript_id ~ "Yes",
    TRUE ~ "No"
  ))

# get probes annotated to canonical transcript promoters
promoter_probes <- annot %>%
  dplyr::filter(Gene_Feature == "promoter" & canonical == "Yes") %>%
  pull(Probe_ID)

# filter methyl data for promoter probes
promoter_methyl <- methyl[methyl$Probe_ID %in% promoter_probes,]

promoter_annot <- annot %>%
  dplyr::filter(Gene_Feature == "promoter" & canonical == "Yes")

promoter_methyl <- promoter_methyl %>%
  left_join(promoter_annot[,c("Probe_ID", "Gene_symbol")])

# calculate mean probe beta value z-scores by plot group
promoter_methyl_df <- promoter_methyl %>%
  gather(key = "Kids_First_Biospecimen_ID_methyl", value = "beta_value", -Probe_ID, -Gene_symbol) %>%
  group_by(Kids_First_Biospecimen_ID_methyl, Gene_symbol) %>%
  summarise(mean_beta_value = mean(beta_value, na.rm = T)) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_methyl,
                                          plot_group)) %>%
  group_by(Gene_symbol, plot_group) %>%
  dplyr::mutate(z_score = scale(mean_beta_value))

# convert to matrix
promoter_zscore_mat <- promoter_methyl_df %>%
  ungroup() %>%
  dplyr::select(-plot_group, -mean_beta_value) %>%
  spread(Kids_First_Biospecimen_ID_methyl, z_score)

# save matrix
saveRDS(promoter_zscore_mat, 
        file.path(results_dir, 
                  "promoter-methyl-zscores.rds"))


# Get all genic probes
gene_probes <- annot %>%
  dplyr::filter(Gene_Feature %in% c("intron", "exon"),
                canonical == "Yes",
                !Probe_ID %in% promoter_probes) %>%
  pull(Probe_ID)

# filter methyl for genic probes
gene_methyl <- methyl[methyl$Probe_ID %in% gene_probes,]

gene_annot <- annot %>%
  dplyr::filter(Gene_Feature %in% c("intron", "exon"),
                canonical == "Yes",
                !Probe_ID %in% promoter_probes) %>%
  distinct(Probe_ID, .keep_all = T)

gene_methyl <- gene_methyl %>%
  left_join(gene_annot[,c("Probe_ID", "Gene_symbol")])

# calculate mean gene body beta value z-scores by plot group
gene_methyl_df <- gene_methyl %>%
  gather(key = "Kids_First_Biospecimen_ID_methyl", value = "beta_value", -Probe_ID, -Gene_symbol) %>%
  group_by(Kids_First_Biospecimen_ID_methyl, Gene_symbol) %>%
  summarise(mean_beta_value = mean(beta_value, na.rm = T)) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID_methyl,
                                   plot_group)) %>%
  group_by(Gene_symbol, plot_group) %>%
  dplyr::mutate(z_score = scale(mean_beta_value))

# convert to matrix
gene_zscore_mat <- gene_methyl_df %>%
  ungroup() %>%
  dplyr::select(-plot_group, -mean_beta_value) %>%
  spread(Kids_First_Biospecimen_ID_methyl, z_score)

saveRDS(gene_zscore_mat, 
        file.path(results_dir, "gene-methyl-zscores.rds"))

# extract H3-wt HGG ids and obtain beta values
hgg_ids <- hist %>%
  dplyr::filter(plot_group == "Other high-grade glioma",
                grepl("H3 wildtype", molecular_subtype)) %>%
  pull(Kids_First_Biospecimen_ID_methyl)

hgg_methyl <- methyl[,colnames(methyl) %in% c("Probe_ID", hgg_ids)]

# Save H3-wt HGG beta value matrix
write_tsv(hgg_methyl, 
        file.path(results_dir, 
                  "hgg-methyl-beta-value.tsv.gz"))


# CPG probe beta values

plp <- read_tsv(plp_file) %>%
  bind_rows(read_tsv(plp_lowVAF_file))

plp_sv <- read_tsv(plp_sv_file)

# create matrix of CPG-annotated probe beta values
cpg_methyl <- methyl %>%
  dplyr::filter(Probe_ID %in% annot$Probe_ID) %>%
  left_join(annot[,c("Probe_ID", "Gene_symbol", "Gene_Feature")]) %>%
  dplyr::filter(Gene_symbol %in% c(plp$gene_symbol_vep, plp_sv$Hugo_Symbol_cpg)) %>%
  gather(key = "Kids_First_Biospecimen_ID_methyl", value = "beta_value", -Probe_ID, -Gene_symbol, -Gene_Feature)

saveRDS(cpg_methyl, 
        file.path(results_dir,
                  "cpg-methyl-beta-values.rds"))

# save annotation df with canonical status
write_tsv(annot, 
          file.path(results_dir,
                    "annot-with-canonical.tsv"))
