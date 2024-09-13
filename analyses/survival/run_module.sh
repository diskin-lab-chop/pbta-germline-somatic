#!/bin/bash

set -e
set -o pipefail

# Format survival data
Rscript -e "rmarkdown::render('01-format-survival.Rmd')"

# Run survival models within histologies by CPG P/LP status
Rscript -e "rmarkdown::render('02-run-survival-plp-status.Rmd')"

# Plot survival models
Rscript --vanilla 03-plot-survival-plp-status.R

# Assess HGG survival by DNA repair gene P/LP status
Rscript -e "rmarkdown::render('04-run-survival-dna-repair.Rmd')"

# Assess braf fusion breakpoint groups and survival by P/LP carrier status
Rscript -e "rmarkdown::render('05-braf-fusion-breakpoint-dist-survival.Rmd')"

# Generate survival summary stats plot
Rscript --vanila 06-survival-summary.R

# Assess P-LP carrier distribution across MB molecular and methylation subtypes
Rscript -e "rmarkdown::render('07-mb-plp-distribution.Rmd')"