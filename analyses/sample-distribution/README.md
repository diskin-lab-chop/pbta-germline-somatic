# Plot CPG P-LP distribution across cohort and tumor histologies

This module creates distribution plots of CPG P-LP variants across the PBTA germline cohort and within tumor histologies (i.e., plot groups). Additional code is included to create plots with PBTA X01 samples for R03 grant. 

## Usage

`bash run_module.sh` 

## Folder content 

1. `01-create-cancer-distribution-plot.R` runs plotting script

## Directory Structure
```
.
├── 01-create-cancer-distribution-plot.R
├── README.md
├── input
│   └── 838_OpenPBTAcohort_1720_X01cohort.tsv
├── plots
│   ├── CPG-PLP-freq-by-histology.png
│   ├── histology-CPG-PLP-count.tiff
│   ├── histology-CPG-PLP-freq.tiff
│   ├── histology-distribution.tiff
│   ├── plp_all_genes-calls.tiff
│   ├── plp_cpgs-calls.tiff
│   └── x01-histology-distribution.tiff
└── run_module.sh
```
  
