# Generate histologies file for OpenPBTA germline project patients

This module creates a cohort-specific histologies file including participant IDs, normal and tumor WGS/WXS sample IDs, and pathology, clinical, and demographic data. 

## Usage

`bash run_module.sh`

## Folder content 

1. `01-collapse-tumor-histologies.R` reads in OpenPedCan V11 histologies file, subsets for patients in germline cohort, and adds additional metadata. 

2. `02-cohort-summary-plots.R` generates cohort summary plots 

## Directory structure
```
.
├── 01-collapse-tumor-histologies.R
├── 02-cohort-summary-plots.R
├── README.md
├── input
│   ├── DEI_CBTN-PNOC_rerun.somalier-ancestry.tsv
│   ├── condition_NOS_pts.tsv
│   ├── plot-mapping.tsv
│   └── samples_of_interest.txt
├── plots
│   ├── circos-plot.pdf
│   └── demo-piecharts.pdf
├── results
│   ├── germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv
│   ├── germline-primary-plus-tumor-histologies-plot-groups.tsv
│   ├── germline-primary-plus-tumor-histologies.tsv
│   └── plot_group_counts.tsv
└── run_module.sh
```