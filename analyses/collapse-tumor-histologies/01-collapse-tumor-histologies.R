library(tidyverse) 

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "collapse-tumor-histologies")
results_dir <- file.path(analysis_dir, "results")

# If the results directory does not exist, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


v22_hist <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"), guess_max = 3000)

# select all normal BS_ids
germline_ids <- v22_hist %>%
  filter(is.na(pathology_diagnosis), 
         experimental_strategy == "WGS") %>%
  select(Kids_First_Biospecimen_ID) %>%
  unique()

# select all tumors
tumor_ids <- v22_hist %>%
  filter(sample_type == "Tumor",
         composition != "Derived Cell Line",
         experimental_strategy == "WGS") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID,
                sample_id_tumor =  sample_id)

primary_tumors <- v22_hist %>%
  filter(sample_type == "Tumor",
         composition != "Derived Cell Line",
         experimental_strategy == "WGS",
         tumor_descriptor == "Initial CNS Tumor") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID,
                sample_id_tumor =  sample_id)

germline_ids_meta <- v22_hist %>%
  filter(Kids_First_Biospecimen_ID %in% germline_ids$Kids_First_Biospecimen_ID) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID,
                sample_id_normal =  sample_id) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal)

combined <- germline_ids_meta %>%
  left_join(tumor_ids) %>%
  group_by(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal) %>%
  summarise(Kids_First_Biospecimen_ID_tumor = str_c(unique(Kids_First_Biospecimen_ID_tumor), collapse = "; "),
            tumor_descriptor = str_c(unique(tumor_descriptor), collapse = "; "),
            cancer_group = str_c(unique(cancer_group), collapse = "; "),
            broad_histology = str_c(unique(broad_histology), collapse = "; "),
            molecular_subtype = str_c(unique(molecular_subtype), collapse = "; "),) %>%
  write_tsv(file.path(results_dir, "tumor-histologies-collapsed-by-germline.tsv"))

combined_primary_only <- germline_ids_meta %>%
  left_join(primary_tumors) %>%
  group_by(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal) %>%
  summarise(Kids_First_Biospecimen_ID_tumor = str_c(unique(Kids_First_Biospecimen_ID_tumor), collapse = "; "),
            tumor_descriptor = str_c(unique(tumor_descriptor), collapse = "; "),
            cancer_group = str_c(unique(cancer_group), collapse = "; "),
            broad_histology = str_c(unique(broad_histology), collapse = "; "),
            molecular_subtype = str_c(unique(molecular_subtype), collapse = "; "),) %>%
  write_tsv(file.path(results_dir, "primary-tumor-histologies-collapsed-by-germline.tsv"))

  
