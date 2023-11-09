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

# load v12 OpenPedCan histologies file + additional samples
v12_hist <- read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000) %>%
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
                                     # neuronal glial NOS to mixed neuronal-glial tumor
                                     pathology_diagnosis == "Glial-neuronal tumor NOS" ~ "Neuronal and mixed neuronal-glial tumor",
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
BS_AFBPM6CN <- v12_hist %>%
  filter(Kids_First_Biospecimen_ID == "BS_AFBPM6CN") %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

# select all normal BS_ids of interest
germline_ids <- read_lines(file.path(input_dir, "samples_of_interest.txt"))

# add independent specimen file to get primary plus tumor bs ids
tumor_ids <- read_tsv(file.path(data_dir, "independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv")) %>%
  bind_rows(BS_AFBPM6CN) %>%
  select(Kids_First_Biospecimen_ID) %>%
  inner_join(v12_hist[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "sample_id", "tumor_descriptor", "pathology_diagnosis", "pathology_free_text_diagnosis", "cancer_group", "broad_histology", "molecular_subtype")]) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID,
                sample_id_tumor = sample_id)

# gather additional germline metadata
germline_ids_meta <- v12_hist %>%
  filter(Kids_First_Biospecimen_ID %in% germline_ids) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID,
                sample_id_normal =  sample_id) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sample_id_normal, germline_sex_estimate) %>%
  distinct()

path_dx <- v12_hist %>%
  filter(!is.na(pathology_diagnosis),
         cohort == "PBTA") %>%
  select(Kids_First_Biospecimen_ID, pathology_diagnosis, pathology_free_text_diagnosis) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID) %>%
  unique() 

# combine germline + tumor ids
combined <- germline_ids_meta %>%
  left_join(tumor_ids) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies.tsv"))

# Read in HOPE project histologies file for updated histology/subtyping info for select patients
hope_hist <- read_tsv(file.path(input_dir, "Hope-GBM-histologies.tsv")) %>%
  dplyr::filter(sample_type == "Tumor") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID)

# add in hope broad_histology, cancer_group, and molecular subtype for HOPE pts
combined <- combined %>%
  left_join(hope_hist[,c("Kids_First_Biospecimen_ID_tumor", "broad_histology",
                         "cancer_group", "molecular_subtype")],
            by = "Kids_First_Biospecimen_ID_tumor",
            suffix = c("_opc", "_hope")) %>%
  # Take HOPE calls when available, OPC otherwise
  dplyr::mutate(broad_histology = case_when(
    !is.na(broad_histology_hope) ~ broad_histology_hope,
    TRUE ~ broad_histology_opc
    ),
    cancer_group = case_when(
      !is.na(cancer_group_hope) ~ cancer_group_hope,
      TRUE ~ cancer_group_opc
    ),
    molecular_subtype = case_when(
      !is.na(molecular_subtype_hope) ~ molecular_subtype_hope,
      TRUE ~ molecular_subtype_opc
    )) %>%
  # change broad histology to DAOT for HGG pts, NMNGT for GNT pt
  dplyr::mutate(broad_histology = case_when(
    Kids_First_Participant_ID %in% c("PT_T8V9ES93", "PT_P5HHJJPH", "PT_3AWKWXEV") ~ "Diffuse astrocytic and oligodendroglial tumor",
    Kids_First_Participant_ID == "PT_SYHB12RN" ~ "Neuronal and mixed neuronal-glial tumor",
    TRUE ~ broad_histology
  ))

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
tumor_clin_meta <- v12_hist %>%
  filter(Kids_First_Biospecimen_ID %in% tumor_ids$Kids_First_Biospecimen_ID_tumor) %>%
  select(sample_id, tumor_descriptor, race, ethnicity, cancer_predispositions, age_at_diagnosis_days, age_last_update_days, OS_days, OS_status, 
         extent_of_tumor_resection, CNS_region, molecular_subtype) %>%
  dplyr::rename(sample_id_tumor = sample_id)

# Most recent event-free survival (EFS) data must be pulled from base histologies file
hist_base <- read_tsv(file.path(data_dir, "histologies-base.tsv"),
                      guess_max = 100000) %>%
  distinct(Kids_First_Participant_ID, .keep_all = T)

# Add ancestry prediction results from Somalier
ancestry <- read_tsv(file.path(input_dir, "DEI_CBTN-PNOC_rerun.somalier-ancestry.tsv")) %>%
  rename(Kids_First_Biospecimen_ID_normal = `#sample_id`)

# Create final histologies file and write to output
final_hist <- combined_map %>%
  left_join(tumor_clin_meta) %>%
  left_join(hist_base[,c("Kids_First_Participant_ID", "EFS_days", "EFS_event_type")]) %>%
  left_join(ancestry[,c("Kids_First_Biospecimen_ID_normal", "predicted_ancestry", "PC1", "PC2", "PC3", "PC4", "PC5")])

predispositions <- read_tsv(file.path(input_dir, "2023-10-31-histologies-base-predispositions.tsv")) %>%
  filter(sample_type == "Tumor") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID)

final_hist <- final_hist %>%
  left_join(predispositions[,c("Kids_First_Biospecimen_ID_tumor", "cancer_predispositions")],
            by = "Kids_First_Biospecimen_ID_tumor",
            suffix = c("_old", "_new")) %>%
  dplyr::mutate(cancer_predispositions = coalesce(cancer_predispositions_new,
                                                  cancer_predispositions_old)) %>%
  dplyr::select(-cancer_predispositions_new, -cancer_predispositions_old) %>%
  write_tsv(file.path(results_dir, "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv"))



