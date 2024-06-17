# Generate demographic and clinical summary statistics among germline P/LP variant carriers in PBTA germline cohort

This module assesses demographic and clinical summary stats among P/LP carriers in PBTA germline cohort and within tumor histologies 

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
├── input
├── plots
│   ├── plot_group_plp_carrier_enr_heatmap.pdf
│   └── predicted-ancestry-pca.pdf
├── results
│   ├── cohort-summary-table-manuscript.tsv
│   ├── demo-clin-pvals-all.tsv
│   ├── demo-clin-pvalues-by-histology.xlsx
│   ├── demo-clin-stats-all.tsv
│   └── demo-clin-stats-by-histology.xlsx
├── run_module.sh
└── util
    ├── heatmap_function.R
    └── summary_functions.R
```