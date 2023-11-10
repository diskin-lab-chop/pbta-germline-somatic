#!/bin/bash

set -e
set -o pipefail

# run incidental findings analysis
Rscript -e "rmarkdown::render('01-incidental-findings.Rmd')"

# plot predisposition findings
Rscript -e "rmarkdown::render('02-plot-predispositions.Rmd')"