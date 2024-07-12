# Characterization of germline P-LP variation in DNA repair genes 

This module summarizes germline pathogenic variation in known DNA repair genes, and further compares mutational signatures and pathway expression in HGG tumors with vs. without P-LP variants in DNA repair pathway genes. 

## Usage

`bash run_dna_repair_variant_analysis.sh` 

## Folder content 

1. `01-format-survival.Rmd` Format molecular subgroups and survival metrics for model generation
2. `02-run-survival-plp-status.Rmd` generate survival models by P/LP carrier status
3. `03-plot-survival-plp-status.R` plot survival models
4. `04-run-survival-dna-repair.Rmd` generate and plot survival models by DNA repair P/LP status in pHGGs. 
5. `05-braf-fusion-breakpoint-dist-survival.Rmd` assesses distribution of BRAF fusion breakpoints among P/LP carriers 


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

├── README.md
├── input
├── plots
│   ├── ATRT/
│   ├── BRCA_km_survival_hgg.pdf
│   ├── CPT/
│   ├── CRANIO/
│   ├── EPN/
│   ├── GNG-GNT/
│   ├── HGG/
│   ├── LGG/
│   ├── MB/
│   ├── MES/
│   ├── MMR_km_survival_hgg.pdf
│   ├── MNG/
│   ├── NFP/
│   ├── forest_EFS_lgg_braf_fusion_1509_resection_plp.pdf
│   ├── forest_EFS_lgg_braf_fusion_1609_resection_plp.pdf
│   ├── forest_EFS_lgg_braf_fusion_resection_group.pdf
│   ├── forest_hgg_add_EFS_subtype_age_repairPLP.pdf
│   ├── forest_hgg_add_OS_subtype_age_repairPLP.pdf
│   ├── forest_hgg_int_EFS_subtype_age_repairPLP.pdf
│   ├── forest_hgg_int_OS_subtype_age_repairPLP.pdf
│   ├── km_EFS_lgg_braf_fusion_1509_plp.pdf
│   ├── km_EFS_lgg_braf_fusion_1609_plp.pdf
│   ├── km_EFS_lgg_braf_fusion_plp.pdf
│   ├── plp_carrier_diagnosis_enr_heatmap.pdf
│   ├── plp_carrier_fusion_enr_heatmap.pdf
│   ├── plp_carrier_region_enr_heatmap.pdf
│   └── plp_carrier_resection_enr_heatmap.pdf
├── results
│   ├── ATRT/
│   ├── CPT/
│   ├── CRANIO/
│   ├── DMG/
│   ├── EPN/
│   ├── GNG-GNT/
│   ├── HGG/
│   ├── LGG/
│   ├── MB/
│   ├── MES/
│   ├── MNG/
│   ├── NFP/
│   ├── coxph_add_EFS_lgg_braf_fusion_1509_resection_plp.RDS
│   ├── coxph_add_EFS_lgg_braf_fusion_1609_resection_plp.RDS
│   ├── coxph_add_EFS_lgg_braf_fusion_resection_group.RDS
│   ├── coxph_add_hgg_EFS_subtype_age_repairPLP.RDS
│   ├── coxph_add_hgg_OS_subtype_age_repairPLP.RDS
│   ├── coxph_int_hgg_EFS_subtype_age_repairPLP.RDS
│   ├── coxph_int_hgg_OS_subtype_age_repairPLP.RDS
│   ├── germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv
│   └── subtypes-for-survival.tsv
├── run_module.sh
└── util
    └── survival_models.R
```
