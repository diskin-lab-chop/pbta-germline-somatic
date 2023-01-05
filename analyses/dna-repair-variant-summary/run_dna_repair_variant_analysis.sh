#!/bin/bash

set -e
set -o pipefail

# Identify samples with germline variants in DNA repair genes
Rscript -e "rmarkdown::render('01-get-dna-repair-variant-samples.Rmd')"

# Compare mutational signature exposures between samples with/without germline DNA repair variants
Rscript -e "rmarkdown::render('02-mutational-signatures.Rmd')"

# Compare GSVA scores between samples with/without germline DNA repair variants
Rscript -e "rmarkdown::render('03-gsva.Rmd')"

# Identify samples with germline variants and somatic mutations in DNA repair pathways
Rscript -e "rmarkdown::render('04-two-hits.Rmd')"
