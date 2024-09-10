#!/bin/bash

set -e
set -o pipefail

# Format methylation data
Rscript --vanilla 01-prepare-methylation.R

# Assess global differential methylation
R -e "rmarkdown::render('02-global-methylation.Rmd')"
