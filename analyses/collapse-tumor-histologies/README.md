# Generate histologies file for OpenPBTA germline project patients

This module creates a cohort-specific histologies file including participant IDs, normal and tumor WGS/WXS sample IDs, and pathology, clinical, and demographic data. 

## Usage

`Rscript 01-collapse-tumor-histologies.R`

## Folder content 

1. `01-collapse-tumor-histologies.R` reads in OpenPedCan V11 histologies file, subsets for patients in germline cohort, and adds additional metadata. 

2. `input/` directory contains the following files: 
  - `CBTN_PNOC.somalier-ancestry.tsv`; includes predicted major ancestry from somalier runs
  - `plot-mapping.tsv`; defines assignment of broad histologies and cancer groups into `plot_group` used for manuscript figures. 
  - `samples_of_interest.txt`; cohort normal BS IDs. 
  
3. `results/` directory contains the following files: 
  - `germline-primary-plus-tumor-histologies.tsv`; includes all histology data 
  - `germline-primary-plus-tumor-histologies-plot-groups.tsv`; includes all histology data + plot group assignment 
  - `germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv`; includes all histology data + plot group assignment + clinical metadata
  - `plot_group_counts.tsv`; reports number of samples in cohort by plot group
  