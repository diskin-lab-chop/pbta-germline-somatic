library(tidyverse) 
library(readxl)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "collapse-tumor-histologies")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")

# If the results directory does not exist, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# load v11 OpenPedCan histologies file
v11_hist <- read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000) %>%
  filter(cohort == "PBTA") %>%
  # update cancer_group to NA if non-tumor
  mutate(cancer_group = case_when(broad_histology == "Non-tumor" ~ NA_character_,
         TRUE ~ as.character(cancer_group)))

# pull cell line info for PT_C2D4JXS1, which has no tumor WGS
BS_AFBPM6CN <- v11_hist %>%
  filter(Kids_First_Biospecimen_ID == "BS_AFBPM6CN") %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

# select all normal BS_ids of interest
germline_ids <- read_lines(file.path(input_dir, "samples_of_interest.txt"))

# add independent specimen file to get primary plus tumor bs ids
tumor_ids <- read_tsv(file.path(data_dir, "independent-specimens.wgs.primary-plus.tsv")) %>%
  bind_rows(BS_AFBPM6CN) %>%
  left_join(v22_hist[,c("Kids_First_Biospecimen_ID", "tumor_descriptor", "cancer_group", "broad_histology", "molecular_subtype")]) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID)
  
# gather additional germline metadata
germline_ids_meta <- v11_hist %>%
  filter(Kids_First_Biospecimen_ID %in% germline_ids) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID,
                sample_id_normal =  sample_id) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal)

path_dx <- v11_hist %>%
  filter(!is.na(pathology_diagnosis),
         cohort == "PBTA") %>%
  select(Kids_First_Biospecimen_ID, pathology_diagnosis, pathology_free_text_diagnosis) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID) %>%
  unique()
  
# combine germline + tumor ids
combined <- germline_ids_meta %>%
  left_join(tumor_ids) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies.tsv"))

# add previous CBTN diagnoses
prev <- readxl::read_excel(file.path(input_dir, "PBTA_Germline_801_10102022.xlsx"), sheet = 1) %>%  dplyr::select(Kids_First_Biospecimen_ID, tumor_descriptor, broad_histology, cancer_group, DISEASE_USING, Other_Description_USE) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID,
                tumor_descriptor_RC = tumor_descriptor,
                broad_histology_RC = broad_histology,
                cancer_group_RC = cancer_group)

combined_append <- combined %>%
  left_join(prev) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, Kids_First_Biospecimen_ID_tumor, 
         tumor_descriptor, broad_histology, broad_histology_RC, cancer_group, cancer_group_RC, DISEASE_USING, Other_Description_USE) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies-CBTN-dx.tsv"))

# add cancer/plot group mapping file 
map_file <- read_tsv(file.path(input_dir, "plot-mapping.tsv"))

# add plot mapping file
combined_map <- combined %>%
  left_join(path_dx) %>%
  left_join(map_file, by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal, 
         Kids_First_Biospecimen_ID_tumor, pathology_diagnosis, pathology_free_text_diagnosis, 
         broad_histology, cancer_group, molecular_subtype, plot_group, broad_histology_display,
         broad_histology_hex, cancer_group_abbreviation, cancer_group_hex, broad_histology_order, 
         oncoprint_group, oncoprint_main) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies-plot-groups.tsv"))


# write plot group counts file
plot_groups <- combined_map %>%
  count(plot_group) %>%
  arrange(n) %>%
  write_tsv(file.path(results_dir, "plot_group_counts.tsv"))


