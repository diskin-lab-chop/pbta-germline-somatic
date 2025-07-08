

library(tidyverse)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
table_dir <- file.path(root_dir, "tables")
input_dir <- file.path(table_dir, "input")
output_dir <- file.path(table_dir, "output")

opc_hist_file <- file.path(data_dir, "histologies.tsv")

cohort_hist_file <- file.path(root_dir, "analyses",
                              "survival",
                              "results", 
                              "germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv")

syndrome_snv_file  <- file.path(root_dir, "analyses",
                                "predisposition-variants",
                                "results", 
                                "incidental-findings-predisposition-variants-path-review.tsv")

syndrome_sv_file  <- file.path(root_dir, "analyses",
                                "predisposition-variants",
                                "results", 
                                "incidental-findings-predisposition-structural-variants-path-review.tsv")

plp_file <- file.path(data_dir,
                      "pbta-merged-plp-variants-autogvp-abridged.tsv")

plp_sv_file <- file.path(data_dir,
                      "pbta_germline_svs.tsv")

cpg_file <- file.path(root_dir, "analyses",
                      "oncokb-annotation",
                      "input", "cpg.txt")

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
                `predicted sex` = germline_sex_estimate,
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
                `broad molecular subgroup`, `predicted sex`,
                 race, ethnicity, `cancer predispositions`,
                `age at diagnosis years`, `age last update years`,
                `OS years`, `OS status`, `EFS years`, `EFS status`,
                `EFS event type`, `extent of tumor resection`,
                `CNS region`, `predicted ancestry`)

write.xlsx(cohort_hist,
           file.path(output_dir,
                     "TableS3.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

# Table S5

syndrome_snvs <- read_tsv(syndrome_snv_file) %>%
  dplyr::rename(`Kids First Participant ID` = Kids_First_Participant_ID,
                `Kids First Biospecimen ID normal` = Kids_First_Biospecimen_ID_normal,
                `Tumor Histology` = plot_group,
                `Molecular subtype` = molecular_subtype,
                Chr = chr,
                Start = start,
                `Ref allele` = ref,
                `Alt allele` = alt,
                `Gene symbol` = gene_symbol_vep,
                `AutoGVP call` = autogvp_call,
                `AutoGVP call reason` = autogvp_call_reason,
                `ClinVar variant ID` = clinvar_variant_id,
                `gnomAD non-cancer AF pop max` = gnomad_3_1_1_AF_non_cancer_popmax,
                `Associated syndrome` = associated_syndrome,
                `Clinically reported syndrome` = clinically_reported_syndrome,
                `Incidental finding` = incidental_finding) %>%
  dplyr::select(-sample_id_normal,
                -sample_id_tumor,
                -clinvar_variant_link,
                -path_report_findings)

syndrome_svs <- read_tsv(syndrome_sv_file) %>%
  dplyr::rename(`Kids First Participant ID` = Kids_First_Participant_ID,
                `Kids First Biospecimen ID normal` = Kids_First_Biospecimen_ID_normal,
                `Tumor Histology` = plot_group,
                `Molecular subtype` = molecular_subtype,
                `Gene symbol` = gene_symbol_vep,
                `Associated syndrome` = associated_syndrome,
                `Clinically reported syndrome` = clinically_reported_syndrome,
                `Incidental finding` = incidental_finding) %>%
  dplyr::select(-sample_id_normal,
                -Kids_First_Biospecimen_ID_tumor,
                -sample_id_tumor,
                -path_report_findings)

syndrome_list <- list("Table S5A-SNVs/InDels" = syndrome_snvs,
                      "Table S5B-SVs" = syndrome_svs)

write.xlsx(syndrome_list,
           file.path(output_dir,
                     "TableS5.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

cpgs <- read_lines(cpg_file)

plp <- read_tsv(plp_file) %>%
  dplyr::filter(gene_symbol_vep %in% cpgs) %>%
  dplyr::rename(`Kids First Biospecimen ID normal` = Kids_First_Biospecimen_ID_normal,
                `variant ids` = variant_ids,
                `Gene symbol` = gene_symbol_vep,
                `Variant classification` = variant_classification_vep,
                `Pathogenicity call` = autogvp_call,
                `Pathogenicity call source` = autogvp_call_reason) %>%
  dplyr::select(-clinvar_clinsig_autogvp,
                -clinvar_review_status_autogvp,
                -clinvar_stars,
                -clinvar_variant_id,
                -clinvar_flag,
                -clinvar_variant_link,
                -intervar_evidence,
                -gnomad_3_1_1_AF_non_cancer_popmax,
                -gnomad_3_1_1_FILTER,
                -variant_id)

plp_sv <- read_tsv(plp_sv_file) %>%
  dplyr::rename(`Kids First Biospecimen ID normal` = Kids_First_Biospecimen_ID_normal,
                chr = Chromosome,
                start = Start,
                alt = Type,
                `SV length` = Length,
                `Gene symbol` = Hugo_Symbol_cpg,
                `Pathogenicity call` = Classification,
                `Pathogenicity call source` = Caller) %>%
  dplyr::select(-Hugo_Symbol_all,
                -End)

merged_plp <- bind_rows(plp, plp_sv) %>%
  dplyr::select(`Kids First Biospecimen ID normal`,
                chr,
                start,
                ref,
                alt,
                `SV length`,
                `variant ids`,
                `Gene symbol`,
                `Variant classification`,
                HGVSg, HGVSc, HGVSp,
                `Pathogenicity call`,
                `Pathogenicity call source`)

sega_pts <- cohort_hist %>%
  dplyr::filter(`plot group` == "Low-grade glioma",
                grepl("SEGA", `molecular subtype`)) %>%
  pull(`Kids First Biospecimen ID normal`)

sega_plp <- merged_plp %>%
  dplyr::filter(`Kids First Biospecimen ID normal` %in% sega_pts)

mb_shh_pts <- cohort_hist %>%
  dplyr::filter(`plot group` == "Medulloblastoma",
                grepl("SHH", `molecular subtype`)) %>%
  pull(`Kids First Biospecimen ID normal`)

mb_shh_plp <- merged_plp %>%
  dplyr::filter(`Kids First Biospecimen ID normal` %in% mb_shh_pts)

hgg_h3wt_pts <- cohort_hist %>%
  dplyr::filter(`molecular subtype` == "HGG, H3 wildtype, TP53") %>%
  pull(`Kids First Biospecimen ID normal`)

hgg_h3wt_plp <- merged_plp %>%
  dplyr::filter(`Kids First Biospecimen ID normal` %in% hgg_h3wt_pts)

pb_pts <- cohort_hist %>%
  dplyr::filter(`cancer group` == "Pineoblastoma") %>%
  pull(`Kids First Biospecimen ID normal`)

pb_plp <- merged_plp %>%
  dplyr::filter(`Kids First Biospecimen ID normal` %in% pb_pts)

subtype_list <- list("Table S8A-SEGA" = sega_plp,
                     "Table S8A-MB-SHH" = mb_shh_plp,
                     "Table S8A-HGG-H3-WT-TP53" = hgg_h3wt_plp,
                     "Table S8A-Pineoblastoma" = pb_plp)

write.xlsx(subtype_list,
           file.path(output_dir,
                     "TableS8.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

# Table S9
table_s9_names <- readxl::excel_sheets(file.path(root_dir, 
                                              "analyses",
                                              "demo-clin-stats",
                                              "results",
                                              "demo-clin-stats-by-histology.xlsx"))

# Read all sheets into a named list
table_s9 <- lapply(table_s9_names, function(x) readxl::read_excel(file.path(root_dir, 
                                                               "analyses",
                                                               "demo-clin-stats",
                                                               "results",
                                                               "demo-clin-stats-by-histology.xlsx"), sheet = x))
table_s9_names <- ifelse(table_s9_names == "Mixed neuronal-glial tumor",
                         "Mixed GNT",
                         ifelse(table_s9_names == "Other high-grade glioma",
                                "Other HGG",
                                ifelse(table_s9_names == "Neurofibroma plexiform",
                                       "NFP", table_s9_names)))

letters_vec <- LETTERS[1:length(table_s9_names)]
names(table_s9) <- glue::glue("TableS10{letters_vec}-{table_s9_names}")

write.xlsx(table_s9,
           file.path(output_dir,
                     "TableS9.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)



# Table S10A
cpg_list_gnomad <- read_tsv(file.path(root_dir, 
                                      "analyses",
                                      "cpg-enrichment",
                                      "input",
                                      "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_cpg_pathway_gnomad_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "gnomAD") %>%
  dplyr::select(-population) %>%
  dplyr::select(pathway_id, genes, pathway_name, `control cohort`, everything())

cpg_list_pmbb <- read_tsv(file.path(root_dir, 
                                    "analyses",
                                    "cpg-enrichment",
                                    "input",
                                    "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_cpg_pathway_pmbb_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "PMBB") %>%
  dplyr::select(pathway_id, genes, pathway_name, `control cohort`, everything())

merged_cpg_list <- cpg_list_gnomad %>%
  bind_rows(cpg_list_pmbb)


cpgs_gnomad <- read_tsv(file.path(root_dir, 
                                      "analyses",
                                      "cpg-enrichment",
                                      "input",
                                      "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_gene_gnomad_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "gnomAD") %>%
  dplyr::rename(`Gene symbol` = gene_symbol_vep) %>%
  dplyr::select(-population) %>%
  dplyr::select(`Gene symbol`, `control cohort`, everything())

cpgs_pmbb <- read_tsv(file.path(root_dir, 
                                    "analyses",
                                    "cpg-enrichment",
                                    "input",
                                    "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_gene_pmbb_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "PMBB") %>%
  dplyr::rename(`Gene symbol` = gene_symbol_vep) %>%
  dplyr::select(`Gene symbol`, `control cohort`, everything())

merged_cpgs <- cpgs_gnomad %>%
  bind_rows(cpgs_pmbb)



kegg_gnomad <- read_tsv(file.path(root_dir, 
                                  "analyses",
                                  "cpg-enrichment",
                                  "input",
                                  "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_kegg_msigdb.v2024.1.Hs.symbols_pathway_gnomad_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "gnomAD") %>%
  dplyr::select(-population) %>%
  dplyr::select(pathway_id, genes, pathway_name, `control cohort`, everything())

kegg_pmbb <- read_tsv(file.path(root_dir, 
                                "analyses",
                                "cpg-enrichment",
                                "input",
                                "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_kegg_msigdb.v2024.1.Hs.symbols.txt_pmbb_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "PMBB") %>%
  dplyr::select(pathway_id, genes, pathway_name, `control cohort`, everything())

merged_kegg <- kegg_gnomad %>%
  bind_rows(kegg_pmbb)



repair_gnomad <- read_tsv(file.path(root_dir, 
                                  "analyses",
                                  "cpg-enrichment",
                                  "input",
                                  "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_dna_repair_pathway_gnomad_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "gnomAD") %>%
  dplyr::select(-population) %>%
  dplyr::select(pathway_id, genes, pathway_name, `control cohort`, everything())

repair_pmbb <- read_tsv(file.path(root_dir, 
                                "analyses",
                                "cpg-enrichment",
                                "input",
                                "pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_dna_repair_pathway_pmbb_enrichment.tsv")) %>%
  dplyr::mutate(`control cohort` = "PMBB") %>%
  dplyr::select(pathway_id, genes, pathway_name, `control cohort`, everything())

merged_repair <- repair_gnomad %>%
  bind_rows(repair_pmbb)


cpg_enr_by_hist <- read_tsv(file.path(root_dir, 
                                      "analyses",
                                      "cpg-enrichment",
                                      "results",
                                      "hist-cpg-plp-enr-pbta-vs-pmbb-gnomad.tsv")) %>%
  dplyr::rename(n_plp_group = n_plp,
                n_no_plp_group = n_no_plp,
                `control cohort` = cohort) %>%
  dplyr::select(pathway_name, plot_group, `control cohort`, everything()) %>%
  dplyr::select(-fraction)


gene_enr_by_hist <- read_tsv(file.path(root_dir, 
                                       "analyses",
                                       "cpg-enrichment",
                                       "results",
                                       "hist-gene-plp-enr-pbta-vs-pmbb-gnomad.tsv")) %>%
  dplyr::rename(n_plp_group = n_plp,
                n_no_plp_group = n_no_plp,
                `control cohort` = cohort,
                `Gene symbol` = gene_symbol_vep) %>%
  dplyr::select(`Gene symbol`, plot_group, `control cohort`, everything()) %>%
  dplyr::select(-fraction,
                -hist_gene)


repair_enr_by_hist <- read_tsv(file.path(root_dir, 
                                         "analyses",
                                         "cpg-enrichment",
                                         "results",
                                         "Knijnenburg_repair_pathways-enrichment-by-hist-PBTA-vs-control.tsv")) %>%
  dplyr::rename(n_plp_group = n_plp,
                n_no_plp_group = n_no_plp,
                `control cohort` = cohort) %>%
  dplyr::select(pathway_name, plot_group, `control cohort`, everything()) %>%
  dplyr::select(-fraction,
                -hist_pathway)


enr_list <- list("TableS10-All CPGs" = merged_cpg_list,
                 "TableS10-CPGs" = merged_cpgs,
                 "TableS10-KEGG pathways" = merged_kegg,
                 "TableS10-DNA repair pathways" = merged_repair,
                 "TableS10-All CPGs by hist" = cpg_enr_by_hist,
                 "TableS10-CPGs by hist" = gene_enr_by_hist,
                 "TableS10-DNA repair by hist" = repair_enr_by_hist)

write.xlsx(enr_list,
           file.path(output_dir,
                     "TableS10.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)


# Table S11
merged_snv_alterations <- read_tsv(file.path(root_dir, "analyses",
                                             "two-hits",
                                             "results",
                                             "plp_snv_indel_somatic_alterations_merged.tsv")) %>%
  dplyr::rename(`Gene symbol` = Hugo_Symbol,
                `VAF tumor` = VAF,
                `Variant classification tumor` = Variant_Classification,
                `TPM z-score` = expr_zscore,
                `Proteomics z-score` = expr_zscore_proteo)

merged_sv_alterations <- read_tsv(file.path(root_dir, "analyses",
                                            "two-hits",
                                            "results",
                                            "plp_sv_somatic_alterations_merged.tsv")) %>%
  dplyr::rename(`Gene symbol` = Hugo_Symbol,
                `VAF tumor` = VAF,
                `Variant classification tumor` = Variant_Classification,
                `TPM z-score` = expr_zscore,
                `Proteomics z-score` = expr_zscore_proteo)  

merged_alterations_list <- list("TableS11A-SNV/InDels" = merged_snv_alterations,
                                "TableS11B-SNV/InDels" = merged_sv_alterations)

write.xlsx(merged_alterations_list,
           file.path(output_dir,
                     "TableS11.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

# Table S12

methyl_df <- read_tsv(file.path(root_dir, "analyses",
                                "methylation",
                                "results",
                                "pbta-germline-methylation-beta-zscores-plus-expr.tsv")) %>%
  dplyr::filter(abs(sample_beta_zscore_group) > 2) %>%
  dplyr::rename(`TPM z-score` = expr_zscore)

write.xlsx(methyl_df,
           file.path(output_dir,
                     "TableS12.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

# Table S13

splice_df <- read_tsv(file.path(root_dir,
                                "analyses",
                                "alternative-splicing",
                                "results",
                                "plp-variant-proximal-splicing-event-psi.tsv")) %>%
  dplyr::rename(`splice event region` = region,
                subgroup = plot_group,
                `include junction counts` = IJC_sample,
                `skip junction counts` = SJC_sample,
                `Percent spliced in (PSI)` = PSI_sample,
                `PSI diff` = PSI_diff_group,
                `PSI, other subgroup` = PSI_other_group,
                `minimum variant-splice event distance (bp)` = min_dist
                
  ) %>%
  dplyr::select(Hugo_Symbol, Kids_First_Biospecimen_ID_rna,
                subgroup, splicing_case, chr,
                `splice event region`, region_start, 
                region_end, upstream_exon_start,
                upstream_exon_end, downstream_exon_start,
                downstream_exon_end, event_in_gtf,
                `include junction counts`, 
                `skip junction counts`,
                `Percent spliced in (PSI)`,
                `PSI, other subgroup`,
                `PSI diff`, sample_PSI_zscore,
                plp_variant_position, ref, alt, 
                variant_classification_vep,
                clinvar_variant_id,
                strand,
                `minimum variant-splice event distance (bp)`)

write.xlsx(splice_df,
           file.path(output_dir,
                     "TableS13.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

# Table S14

hgg_mutsigs_diff <- read_tsv(file.path(root_dir,
                                       "analyses",
                                       "dna-repair-variant-summary",
                                       "results",
                                       "HGG, H3 wildtype",
                                       "cosmicv3_sig_diff_MMR_vs_other.tsv"))


sig_map <- read_tsv(file.path(input_dir, 
                              "sbs_v3_map.tsv"))

table_s14 <- sig_map %>%
  left_join(hgg_mutsigs_diff,
            by = c("Signature" = "signature")) %>%
  arrange(pvalue)

write.xlsx(table_s14,
           file.path(output_dir,
                     "TableS14.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

tmb_methyl_corr <- read_tsv(file.path(root_dir, 
                                      "analyses",
                                      "methylation",
                                      "results",
                                      "tmb-methylation-correlation-by-hist.tsv")) %>%
  dplyr::rename(`TMB-mean methylation Pearson r` = pearson_r,
                `TMB-mean methylation Pearson p-value` = pearson_p) %>% 
  arrange(`TMB-mean methylation Pearson p-value`)

write.xlsx(tmb_methyl_corr,
           file.path(output_dir,
                     "TableS15.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)

# Table S16

surv_df <- read_tsv(file.path(root_dir, 
                              "analyses",
                              "survival",
                              "results",
                              "median-survival-by-ancestry-cancer-group.tsv")) %>%
  dplyr::arrange(type, group) %>%
  dplyr::select(-group_plus_n)

write.xlsx(surv_df,
           file.path(output_dir,
                     "TableS16.xlsx"),
           overwrite = TRUE,
           keepNA = TRUE,
           rowNames = FALSE)


