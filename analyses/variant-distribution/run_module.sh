#!/bin/bash

set -e
set -o pipefail

# Run variant distribution 
Rscript -e "rmarkdown::render('01-create-variant-distribution-plots.Rmd')"
