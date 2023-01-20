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

# load v11_plus OpenPedCan histologies file + additional samples
v11_hist <- read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000) %>%
  filter(cohort == "PBTA") %>%
  # update choroid plexus from benign to tumor category
  mutate(broad_histology = case_when(pathology_diagnosis == "Choroid plexus papilloma" ~ "Choroid plexus tumor",
                                     pathology_diagnosis == "Choroid plexus carcinoma" ~ "Choroid plexus tumor",
                                     pathology_free_text_diagnosis == "choroid plexus cyst" ~ "Choroid plexus tumor",
                                     # make cavernoma benign tumor
                                     pathology_diagnosis == "Cavernoma" ~ "Benign tumor",
                                     # meningiomatosis --> meningioma
                                     pathology_free_text_diagnosis == "meningioangiomatosis" ~ "Meningioma",
                                     # MPNST --> NF
                                     grepl("MPNST", pathology_diagnosis) ~ "Tumor of cranial and paraspinal nerves",
                                     # move these from benign to other tumor
                                     pathology_free_text_diagnosis %in% c("fibrous dysplasia",
                                                                          "osteoblastoma",
                                                                          "pituitary macroadenoma",
                                                                          "pituitary adenoma",
                                                                          "prolactinoma",
                                                                          "perineuroma") ~ "Other tumor",
                                     broad_histology == "Other" ~ "Other tumor",
                                     # pineoblastoma
                                     pathology_free_text_diagnosis == "pnet - pineoblastoma with calcification" ~ "Tumor of pineal region",
                                     # move oligo II to LGG
                                     pathology_free_text_diagnosis == "oligodendroglioma who ii" ~ "Low-grade astrocytic tumor",
                                     pathology_diagnosis == "Oligodendroglioma" ~ "Low-grade astrocytic tumor",
                                     TRUE ~ as.character(broad_histology)),
         # update cancer_group to NA if non-tumor
         cancer_group = case_when(broad_histology %in% c("Non-tumor", "Benign Tumor") ~ NA_character_,
                                  grepl("MPNST", pathology_diagnosis) ~ "Neurofibroma/Plexiform",
                                  pathology_free_text_diagnosis == "meningioangiomatosis" ~ "Meningioma",
                                  # move these from benign to other tumor
                                  pathology_free_text_diagnosis %in% c("fibrous dysplasia",
                                                                       "osteoblastoma",
                                                                       "pituitary macroadenoma",
                                                                       "pituitary adenoma",
                                                                       "prolactinoma",
                                                                       "perineuroma") ~ "Other tumor",
                                  broad_histology == "Other tumor" ~ "Other tumor",
                                  # pineoblastoma
                                  pathology_free_text_diagnosis == "pnet - pineoblastoma with calcification" ~ "Pineoblastoma",
                                  TRUE ~ as.character(cancer_group)),

  )
                                     
  

# pull cell line info for PT_C2D4JXS1, which has no tumor WGS
BS_AFBPM6CN <- v11_hist %>%
  filter(Kids_First_Biospecimen_ID == "BS_AFBPM6CN") %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

# select all normal BS_ids of interest
germline_ids <- read_lines(file.path(input_dir, "samples_of_interest.txt"))

# add independent specimen file to get primary plus tumor bs ids
tumor_ids <- read_tsv(file.path(data_dir, "independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv")) %>%
  bind_rows(BS_AFBPM6CN) %>%
  select(Kids_First_Biospecimen_ID) %>%
  inner_join(v11_hist[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "tumor_descriptor", "pathology_diagnosis", "pathology_free_text_diagnosis", "cancer_group", "broad_histology", "molecular_subtype")]) %>%
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

# add plot mapping file and old plot groups
combined_map <- combined %>%
  left_join(path_dx) %>%
  left_join(prev[,c("Kids_First_Biospecimen_ID_normal", "DISEASE_USING","Other_Description_USE")]) %>%
  left_join(map_file, by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal, 
         Kids_First_Biospecimen_ID_tumor, pathology_diagnosis, pathology_free_text_diagnosis, 
         broad_histology, cancer_group, plot_group, Other_Description_USE, DISEASE_USING, molecular_subtype, broad_histology_display,
         broad_histology_hex, cancer_group_abbreviation, cancer_group_hex, broad_histology_order, 
         oncoprint_group, oncoprint_main) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies-plot-groups.tsv"))


# write plot group counts file
plot_groups <- combined_map %>%
  count(plot_group) %>%
  arrange(n) %>%
  write_tsv(file.path(results_dir, "plot_group_counts.tsv"))

# previous plot groups
plot_groups_previous <- combined_map %>%
  count(Other_Description_USE) %>%
  arrange(n) %>%
  write_tsv(file.path(results_dir, "previous_plot_group_counts.tsv"))
