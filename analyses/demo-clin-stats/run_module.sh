#!/bin/bash

set -e
set -o pipefail

# Get demo-clin stats by CPG P/LP status
Rscript -e "rmarkdown::render('01-demo-clin-stats.Rmd')"

# Plot somalier predicted ancestry
Rscript --vanilla 02-plot-ancestry.R