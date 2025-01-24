#!/bin/bash

set -e
set -o pipefail

# Run AS event z-score calculation script
Rscript -e "rmarkdown::render('01-calculate-PSI-zscores.Rmd')"

# Run PSI-TPM correlation assessment script
Rscript -e "rmarkdown::render('02-splicing-expression.Rmd')"

# Identify P/LP-variant associated alternative splicing events
Rscript -e "rmarkdown::render('03-plp-proximal-alternative-splicing.Rmd')"