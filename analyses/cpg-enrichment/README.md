# Assess germline PLP variant enrichment in OpenPBTA germline cohort

This module calculates enrichment of P-LP variants in cancer predisposition genes within PBTA germline cohort relative to two cancer-free control cohorts (Penn Med Biobank and gnomAD), and plots significant results. 

## Usage

`Rscript -e "rmarkdown::render('01-run-enrichment.Rmd', clean = TRUE)"`
`Rscript 02-plot-enrichment.R`

## Folder content 

1. `01-run-enrichment.Rmd` runs Fishers exact tests on all cancer predisposition genes to calculate enrichment of P-LP variants in PBTA germline cohort relative to PMBB and gnomAD. Enrichment is also calculated for genes within tumor histologies (i.e., plot groups). 

2. `02-plot-enrichment.R` plots significant enrichment results

3. `input/` directory contains the following files: 
  - `PBTA-CPG-PLP-enrichment-PMBB.tsv`; contains counts of P-LP variants in the PMBB cohort for each CPG
  - `PBTA-CPG-PLP-enrichment-gnomAD.tsv`; contains counts of P-LP variants in the gnomAD cohort for each CPG

4. `results/` directory contains the following files: 
  - `PBTA-cpg-plp-enrichment-PMBB.tsv`; enrichment in PBTA vs. PMBB 
  - `PBTA-cpg-plp-enrichment-by-plotGroup-PMBB.tsv`; enrichment in PBTA plot groups vs. PMBB 
  - `PBTA-cpg-plp-enrichment-by-plotGroup-gnomAD.tsv`; enrichment in PBTA plot groups vs. gnomAD 
  - `PBTA-cpg-plp-enrichment-gnomAD.tsv`; enrichment in PBTA vs. gnomAD 
  - `cpg-plp-enrichment-merged.tsv`; merged enrichment results (all cohort)
  - `plot-group-cpg-plp-enrichment-merged.tsv`; merged enrichment results 
  
5. `plots/` directory contains the following files: 
  - `sig-CPG-enrichment-PBTA-vs-control.tiff`; plots of significant enrichments across cohort
  - `sig-hist-CPG-enrichment-PBTA-vs-control.png` plots of significant enrichments within cohort tumor histologies