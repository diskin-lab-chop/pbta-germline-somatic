#!/bin/bash

set -e
set -o pipefail

# Format survival data
Rscript -e "rmarkdown::render('01-format-survival.Rmd')"

# Run survival models within histologies by CPG P/LP status
Rscript -e "rmarkdown::render('02-run-survival-plp-status.Rmd')"

# Plot survival models
Rscript 03-plot-survival-plp-status.R

# Assess HGG survival by DNA repair gene P/LP status
Rscript -e "rmarkdown::render('04-run-survival-dna-repair.Rmd')"