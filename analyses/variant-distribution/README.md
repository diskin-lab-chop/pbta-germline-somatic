# Plot P-LP variant distribution across cancer predisposition genes

This module creates distribution plots of P-LP variants across 196 cancer predisposition genes (CPGs). Additional code is included to create plots with PBTA X01 samples for R03 grant. 

## Usage

`Rscript -e "rmarkdown::render('01-create-variant-distribution-plots.Rmd')"` 

## Folder content 

1. `01-create-variant-distribution-plots.Rmd` runs plotting script

2. `input/` directory contains the following files:
  - `embryonal-tumor_goi_list.tsv`; list of embryonal tumor driver genes
  - `hgat_goi_list.tsv`; list of hgat driver genes
  - `lgat_goi_list.tsv`; list of lgg tumor driver genes
  - `other_goi_list.tsv`; list of other tumor driver genes

3. `results/` directory contains the following files: 
  - `plp-variants-in-somatic-drivers-not-cpg.tsv`
  
4. `plots/` directory contains the following figures: 
  - Barplots of no. P-LP variants by CPG (n>=2) in 838 PBTA cohort and + X01 cohort (`cpg-variant-distribution-by-histology.tiff`, `cpg-variant-distribution-by-histology_R03.tiff`)
  - Barplot of no. P-LP variants in somatic driver genes that are also CPGs (`somatic_drivers-variant-distribution-by-histology.tiff`)
  - Barplot of no. P-L variants in somatic driver genes that are not CPGs (`somatic_not_cpg-variant-distribution-by-histology.tiff`)

