

library(tidyverse)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
table_dir <- file.path(root_dir, "tables")
output_dir <- file.path(table_dir, "output")

opc_hist_file <- file.path(data_dir, "histologies.tsv")

cohort_hist_file <- file.path(root_dir, "analyses",
                              "survival",
                              "results", 
                              "germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv")

# wrangle data

opc_hist <- read_tsv(opc_hist_file)

cohort_hist <- read_tsv(cohort_hist_file) %>%
  # remove columns to not include in table
  dplyr::select(-contains("PC"),
                -broad_histology_display,
                -broad_histology_hex,
                -cancer_group_abbreviation,
                -plot_group_hex,
                -broad_histology_order,
                -oncoprint_group,
                -cohort_participant_id,
                -cpgPLP_status,
                -OS_days,
                -EFS_days) %>%
  # add sub_cohort
  left_join(opc_hist %>% dplyr::select(Kids_First_Biospecimen_ID,
                                       reported_gender,
                                       sub_cohort),
            by = c("Kids_First_Biospecimen_ID_normal" = "Kids_First_Biospecimen_ID")) %>%
  # convert age at diagnosis and last update to years
  dplyr::mutate(age_at_diagnosis_years = age_at_diagnosis_days/365.25,
                age_last_update_years = age_last_update_days/365.25) %>%
  # rename column headers
  dplyr::rename(`Kids First Participant ID` = Kids_First_Participant_ID,
                `cohort` = sub_cohort,
                `Kids First Biospecimen ID normal` = Kids_First_Biospecimen_ID_normal,
                `sample id normal` = sample_id_normal,
                `Kids First Biospecimen ID tumor` = Kids_First_Biospecimen_ID_tumor,
                `sample id tumor` = sample_id_tumor,
                `tumor event` = tumor_descriptor,
                `pathology diagnosis` = pathology_diagnosis,
                `pathology free text diagnosis` = pathology_free_text_diagnosis,
                `broad histology` = broad_histology,
                `cancer group` = cancer_group,
                `plot group` = plot_group,
                `molecular subtype` = molecular_subtype,
                `broad molecular subgroup` = mol_sub_group,
                `germline sex estimate` = germline_sex_estimate,
                `reported gender` = reported_gender,
                `cancer predispositions` = cancer_predispositions,
                `age at diagnosis years` = age_at_diagnosis_years,
                `age last update years` = age_last_update_years,
                `OS years` = OS_years,
                `OS status` = OS_status,
                `EFS years` = EFS_years,
                `EFS status` = EFS_status,
                `EFS event type` = EFS_event_type,
                `extent of tumor resection` = extent_of_tumor_resection,
                `CNS region` = CNS_region,
                `predicted ancestry` = predicted_ancestry) %>%
  # re-order columns
  dplyr::select(`Kids First Participant ID`, `cohort`, `Kids First Biospecimen ID normal`,
                `sample id normal`, `Kids First Biospecimen ID tumor`,
                `sample id tumor`, `tumor event`, `pathology diagnosis`,
                `pathology free text diagnosis`, `broad histology`,
                `cancer group`, `plot group`, `molecular subtype`,
                `broad molecular subgroup`, `germline sex estimate`,
                `reported gender`, race, ethnicity, `cancer predispositions`,
                `age at diagnosis years`, `age last update years`,
                `OS years`, `OS status`, `EFS years`, `EFS status`,
                `EFS event type`, `extent of tumor resection`,
                `CNS region`, `predicted ancestry`)

# write to output
write_tsv(cohort_hist,
          file.path(output_dir,
                    "Table-S3A.tsv"))

# Create table of patient freq by cohort and plot group
hist_cohort_df <- table(cohort_hist$`plot group`,
                        cohort_hist$cohort)

hist_cohort_df <- hist_cohort_df %>%
  as.data.frame() %>%
  pivot_wider(names_from = Var2,
              values_from = Freq) %>%
  dplyr::rename(`plot group` = Var1)

# write to output
write_tsv(hist_cohort_df,
          file.path(output_dir,
                      "Table-S3B.tsv"))

