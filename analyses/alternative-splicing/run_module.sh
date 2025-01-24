#!/bin/bash

set -e
set -o pipefail

# Run AS event z-score calculation script
Rscript -e "rmarkdown::render('01-calculate-PSI-zscores.Rmd')"

# Run PSI-TPM correlation assessment script
Rscript -e "rmarkdown::render('02-splicing-expression.Rmd')"