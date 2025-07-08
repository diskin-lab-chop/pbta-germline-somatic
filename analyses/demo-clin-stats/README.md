# Generate demographic and clinical summary statistics among germline P/LP variant carriers in PBTA germline cohort

This module assesses demographic and clinical summary stats among P/LP carriers in PBTA germline cohort and within tumor histologies and generates multiple tables for the manuscript.

## Main Table 1
Derived from `demo-clin-stats-cohort.xlsx`

## Supplemental Table S9
Derived from `demo-clin-stats-by-histology.xlsx`


## Usage

`bash run_module.sh` 

## Folder content 

1. `01-demo-clin-stats.Rmd` generate demographic and clinical summary stats

2. `02-plot-ancestry.R` plot Somalier genetic ancestry PCs

## Directory structure
```
.
├── 01-demo-clin-stats.Rmd
├── 01-demo-clin-stats.nb.html
├── 02-plot-ancestry.R
├── README.md
├── plots
│   ├── other_group_plp_carrier_enr.pdf
│   ├── plot_group_plp_carrier_enr_heatmap.pdf
│   ├── predicted-ancestry-pca.pdf
│   └── subtype_plp_carrier_enr.pdf
├── results
│   ├── demo-clin-pvals-all.tsv
│   ├── demo-clin-stats-all.tsv
│   ├── demo-clin-stats-by-histology.xlsx
│   ├── demo-clin-stats-cohort.xlsx
│   └── plp-carrier-subtype-counts.tsv
├── run_module.sh
└── util
    ├── heatmap_function.R
    └── summary_functions.R
```