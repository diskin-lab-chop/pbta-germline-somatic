# Survival analyses in PBTA germline cohort

This module assesses overall and event-free survival in the PBTA germline cohort to determine if germline P-LP variant carriers exhibit different outcomes compared to non-carriers. 

## Usage

`bash run_module.sh` 

## Folder content 

1. `01-format-survival.Rmd` Format molecular subgroups and survival metrics for model generation
2. `02-run-survival-plp-status.Rmd` generate survival models by P/LP carrier status
3. `03-plot-survival-plp-status.R` plot survival models
4. `04-run-survival-dna-repair.Rmd` generate and plot survival models by DNA repair P/LP status in pHGGs. 
5. `05-braf-fusion-breakpoint-dist-survival.Rmd` assesses distribution of BRAF fusion breakpoints among P/LP carriers 
6. `06-survival-summary.R` plots P-LP carrier OS and EFS hazard ratios for each histologya and molecular subtype cohort
7. `07-mb-plp-distribution.Rmd` assessees distribution of P-LP carriers across MB molecular and methylation subtypes
8. `08-mb-survival.Rmd` assesses survival in MB cohort by molecular and methylation subtypes and metastasis status

## Directory structure
```
.
├── 01-format-survival.Rmd
├── 01-format-survival.nb.html
├── 02-run-survival-plp-status.Rmd
├── 02-run-survival-plp-status.nb.html
├── 03-plot-survival-plp-status.R
├── 03-plot-survival-plp-status.html
├── 04-run-survival-dna-repair.Rmd
├── 04-run-survival-dna-repair.nb.html
├── 05-braf-fusion-breakpoint-dist-survival.Rmd
├── 05-braf-fusion-breakpoint-dist-survival.nb.html
├── 06-survival-summary.R
├── 07-mb-plp-distribution.Rmd
├── 07-mb-plp-distribution.nb.html
├── 08-mb-survival.Rmd
├── 08-mb-survival.nb.html
├── README.md
├── input
│   ├── cbtn-all-survival-09122024.txt
│   ├── lgg-braf-fusion-breakpoint-annotation.tsv
│   ├── mb_shh_molecular_subtypes.tsv
│   ├── pbta-lgg-braf-fusion-mixed-location.tsv
│   └── progression-metastasis-data.tsv
├── plots
│   ├── ATRT
│   ├── BRCA_km_survival_hgg.pdf
│   ├── CRANIO
│   ├── DMG
│   ├── EPN
│   ├── GNG-GNT
│   ├── HGG
│   ├── LGG
│   ├── MB
│   ├── MMR_km_survival_hgg.pdf
│   ├── MNG
│   ├── NFP
│   ├── breakpoint_group_region_enr_heatmap.pdf
│   ├── breakpoint_group_resection_enr_heatmap.pdf
│   ├── forest_EFS_lgg_braf_fusion_resection_group.pdf
│   ├── forest_EFS_lgg_braf_fusion_resection_group_plp.pdf
│   ├── forest_EFS_mb_Group3_resection_subtype_metastasis.pdf
│   ├── forest_EFS_mb_Group4_resection_subtype_metastasis.pdf
│   ├── forest_EFS_mb_SHH_resection_subtype_metastasis.pdf
│   ├── forest_EFS_mb_WNT_resection_subtype_metastasis.pdf
│   ├── forest_EFS_mb_g3_resection_methyl_subtype.pdf
│   ├── forest_EFS_mb_g4_resection_methyl_subtype.pdf
│   ├── forest_EFS_mb_resection_subtype_metastasis_plp.pdf
│   ├── forest_OS_mb_Group3_resection_subtype_metastasis.pdf
│   ├── forest_OS_mb_Group4_resection_subtype_metastasis.pdf
│   ├── forest_OS_mb_SHH_resection_subtype_metastasis.pdf
│   ├── forest_OS_mb_g3_resection_methyl_subtype.pdf
│   ├── forest_OS_mb_resection_subtype_metastasis_plp.pdf
│   ├── forest_hgg_add_EFS_subtype_age_repairPLP.pdf
│   ├── forest_hgg_add_OS_subtype_age_repairPLP.pdf
│   ├── forest_hgg_int_EFS_subtype_age_repairPLP.pdf
│   ├── forest_hgg_int_OS_subtype_age_repairPLP.pdf
│   ├── km_EFS_lgg_braf_fusion_group.pdf
│   ├── km_OS_mb_shh_subtype.pdf
│   ├── km_mb_g4_subtype.pdf
│   ├── mb_MB, Group3_plp_carrier_metastasis_enr_heatmap.pdf
│   ├── mb_MB, Group4_plp_carrier_metastasis_enr_heatmap.pdf
│   ├── mb_MB, SHH_plp_carrier_metastasis_enr_heatmap.pdf
│   ├── mb_MB, WNT_plp_carrier_metastasis_enr_heatmap.pdf
│   ├── mb_g3_plp_carrier_subtype_enr_heatmap.pdf
│   ├── mb_g4_plp_carrier_subtype_enr_heatmap.pdf
│   ├── mb_plp_carrier_metastasis_enr_heatmap.pdf
│   ├── mb_shh_plp_carrier_subtype_enr_heatmap.pdf
│   ├── mb_shh_plp_carrier_tp53_enr_heatmap.pdf
│   ├── mb_shh_plp_carrier_type_enr_heatmap.pdf
│   ├── plp_carrier_diagnosis_enr_heatmap.pdf
│   ├── plp_carrier_fusion_enr_heatmap.pdf
│   ├── plp_carrier_mb_subtype_enr_heatmap.pdf
│   └── survival-hr-plp-vs-wt.pdf
├── results
│   ├── ATRT
│   ├── CRANIO
│   ├── DMG
│   ├── EPN
│   ├── GNG-GNT
│   ├── HGG
│   ├── LGG
│   ├── MB
│   ├── MNG
│   ├── NFP
│   ├── coxph_add_EFS_lgg_braf_fusion_resection_group.RDS
│   ├── coxph_add_EFS_lgg_braf_fusion_resection_group_plp.RDS
│   ├── coxph_add_EFS_mb_Group3_resection_metastasis_plp.RDS
│   ├── coxph_add_EFS_mb_Group4_resection_metastasis_plp.RDS
│   ├── coxph_add_EFS_mb_SHH_resection_metastasis_plp.RDS
│   ├── coxph_add_EFS_mb_WNT_resection_metastasis_plp.RDS
│   ├── coxph_add_EFS_mb_g3_resection_methyl_subtype.RDS
│   ├── coxph_add_EFS_mb_g4_methyl_subtype.RDS
│   ├── coxph_add_EFS_mb_resection_subtype_metastasis_plp.RDS
│   ├── coxph_add_OS_mb_Group3_resection_metastasis_plp.RDS
│   ├── coxph_add_OS_mb_Group4_resection_metastasis_plp.RDS
│   ├── coxph_add_OS_mb_SHH_resection_metastasis_plp.RDS
│   ├── coxph_add_OS_mb_g3_resection_methyl_subtype.RDS
│   ├── coxph_add_OS_mb_resection_subtype_metastasis_plp.RDS
│   ├── coxph_add_hgg_EFS_subtype_age_repairPLP.RDS
│   ├── coxph_add_hgg_OS_subtype_age_repairPLP.RDS
│   ├── coxph_int_hgg_EFS_subtype_age_repairPLP.RDS
│   ├── coxph_int_hgg_OS_subtype_age_repairPLP.RDS
│   ├── germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv
│   ├── mb-histologies-plus-shh-methyl-subtypes-metastsis.tsv
│   ├── median-survival-by-ancestry-cancer-group.tsv
│   └── subtypes-for-survival.tsv
├── run_module.sh
└── util
    └── survival_models.R
```
