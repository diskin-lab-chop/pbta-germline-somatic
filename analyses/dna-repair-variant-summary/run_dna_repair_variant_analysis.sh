#!/bin/bash

set -e
set -o pipefail

# Identify samples with germline PLP variants in DNA repair pathway genes
Rscript -e "rmarkdown::render('01-get-dna-repair-variant-samples.Rmd')"

# Compare tumor mutational signatures between samples with and without variants in DNA repair genes
Rscript -e "rmarkdown::render('02-mutational-signatures.Rmd')"

# Compare tumor gsva scores between samples with and without variants in DNA repair genes
Rscript -e "rmarkdown::render('03-run-gsva-comparison.Rmd')"
