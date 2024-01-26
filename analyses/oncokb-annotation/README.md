# Run oncokb annotation on all somatic mutations in germline PBTA cohort corresponding tumors

This module calculates enrichment of P-LP variants in cancer predisposition genes within PBTA germline cohort relative to two cancer-free control cohorts (Penn Med Biobank and gnomAD), and plots significant results. 

## Usage

`ONCOKB=<ONCOKB-key> bash run-oncokb-annotation.sh` 

NOTE: an ONCOKB access token is required when running this bash script 

## Folder content 

1. `run-oncokb-annotation.sh` runs all module scripts in `code/` directory

2. `input/` directory contains the following files:
  - `cpg.txt`; list of 196 cancer predisposition genes

3. `code/` directory contains the following files: 
  - `01-create-subset-maf.Rmd` subset maf file for only 838 germline PBTA patient tumor samples
  - `02-run-oncokb-annotation.sh` runs oncokb annotation. NOTE: a ONCOKB access token is required when running this script. 

4. `results/` directory contains the following files: 
  - `snv-consensus-plus-hotspots-goi.maf.tsv`; maf subsetted for germline PBTA patient tumors
  - `snv-consensus-plus-hotspots-goi-oncokb.maf.tsv`; subsetted maf file with oncokb annotation