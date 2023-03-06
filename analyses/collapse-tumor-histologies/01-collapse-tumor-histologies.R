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
  inner_join(v11_hist[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "sample_id", "tumor_descriptor", "pathology_diagnosis", "pathology_free_text_diagnosis", "cancer_group", "broad_histology", "molecular_subtype")]) %>%
dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID,
              sample_id_tumor = sample_id)

# gather additional germline metadata
germline_ids_meta <- v11_hist %>%
  filter(Kids_First_Biospecimen_ID %in% germline_ids) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID,
                sample_id_normal =  sample_id) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal, germline_sex_estimate) %>%
  distinct()

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


# add cancer/plot group mapping file 
map_file <- read_tsv(file.path(input_dir, "plot-mapping.tsv")) %>%
  distinct()

# add plot mapping file and old plot groups
combined_map <- combined %>%
  left_join(path_dx) %>%
  left_join(map_file, by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal, 
         Kids_First_Biospecimen_ID_tumor, sample_id_tumor, tumor_descriptor, pathology_diagnosis, pathology_free_text_diagnosis, 
         broad_histology, cancer_group, plot_group, molecular_subtype, broad_histology_display,
         broad_histology_hex, cancer_group_abbreviation, plot_group_hex, broad_histology_order, 
         oncoprint_group, germline_sex_estimate) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies-plot-groups.tsv"))

# make sure no duplicate normal ids
combined_map[duplicated(combined_map$Kids_First_Biospecimen_ID_normal),]

# write plot group counts file
plot_groups <- combined_map %>%
  count(plot_group) %>%
  arrange(n) %>%
  write_tsv(file.path(results_dir, "plot_group_counts.tsv"))


# finally, add relevant clinical information to new histologies and plot group file
tumor_clin_meta <- v11_hist %>%
  filter(Kids_First_Biospecimen_ID %in% tumor_ids$Kids_First_Biospecimen_ID_tumor) %>%
  select(sample_id, tumor_descriptor, race, ethnicity, cancer_predispositions, age_at_diagnosis_days, age_last_update_days, PFS_days, OS_days, OS_status, 
         extent_of_tumor_resection, CNS_region, molecular_subtype) %>%
  dplyr::rename(sample_id_tumor = sample_id)
  
ancestry <- read_tsv(file.path(input_dir, "pbta-germline-ancestry-results.tsv")) %>%
  rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID)

final_hist <- combined_map %>%
  left_join(tumor_clin_meta) %>%
  left_join(ancestry[,c("Kids_First_Biospecimen_ID_normal", "predicted_ancestry", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")], by = "Kids_First_Biospecimen_ID_normal") %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv"))





