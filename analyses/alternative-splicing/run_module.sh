#!/bin/bash

set -e
set -o pipefail

# if processed rMATs does not exist, run 00-prepare_rmats.sh
if [ -f "results/splice-events-rmats-cpgs.tsv.gz" ]; then
    echo "Found processed rmats file. Proceeding..."
else
    echo "Processed rmats file does not exist. Running 00-filter-rmats.R..."
    Rscript --vanilla 00-filter-rmats.R
fi

# Run AS event z-score calculation script
Rscript -e "rmarkdown::render('01-calculate-PSI-zscores.Rmd')"

# Run PSI-TPM correlation assessment script
Rscript -e "rmarkdown::render('02-splicing-expression.Rmd')"

# Identify P/LP-variant associated alternative splicing events
Rscript -e "rmarkdown::render('03-plp-proximal-alternative-splicing.Rmd')"