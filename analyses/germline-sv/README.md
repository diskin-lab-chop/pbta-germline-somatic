# Plot P-LP structural variant distribution across cancer predisposition genes

This module creates distribution plots of P-LP structural variants across 214 cancer predisposition genes (CPGs).

## Usage

`Rscript -e "rmarkdown::render('01-sv-summary.Rmd')"` 

## Folder content 

1. `01-sv-summary.Rmd` runs plotting script
  
2. `plots/` directory contains the following figures: 
  - Barplots of no. P-LP structural variants by CPG in 837 PBTA cohort by histology SV type (deletion or duplication) (`sv-plp-distribution-plot.tiff`)

