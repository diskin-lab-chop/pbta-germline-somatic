# Plot CPG P-LP distribution across cohort and tumor histologies

This module creates distribution plots of CPG P-LP variants across the PBTA germline cohort and within tumor histologies (i.e., plot groups). Additional code is included to create plots with PBTA X01 samples for R03 grant. 

## Usage

`Rscript 01-create-cancer-distribution-plot.R` 

## Folder content 

1. `01-create-cancer-distribution-plot.R` runs plotting script

2. `input/` directory contains the following files:
  - `838_OpenPBTAcohort_1720_X01cohort.tsv`; patient metadata for PBTA germline cohort + PBTA X01 cohort

3. `plots/` directory contains the following figures: 
  - Lollipop plots of CPG P-LP frequency by plot group in 837 cohort and + X01 (`CPG-PLP-freq-by-histology.png`, `CPG-PLP-freq-by-histology_r03.png`)
  - Barplots of no. and perc. patients with CPG P-LP variant by plot group (`histology-CPG-PLP-count.tiff`, `histology-CPG-PLP-freq.tiff`)
  - Barplots of no. patients by plot group in 837 PBTA germline cohort and + X01 cohort (`histology-distribution.tiff`, `x01-histology-distribution.tiff`)
  - Lollipop plot of number of P-LP variant calls in all genes and CPGs, in 837 cohort and + X01 cohort (`plp_all_genes-calls.tiff`, `plp_cpgs-calls.tiff`, `plp_cpgs_r03-calls.tiff`)
  
  
