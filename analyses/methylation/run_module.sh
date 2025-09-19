#!/bin/bash

set -e
set -o pipefail

# Format methylation data
Rscript --vanilla 01-prepare-methylation.R

# Assess global differential methylation
R -e "rmarkdown::render('02-global-methylation.Rmd')"

# Assess probe differential methylation
R -e "rmarkdown::render('03-calculate-probe-zscores.Rmd')"

# Plot HGG H3 wt promoter methylation
R -e "rmarkdown::render('04-plot-hgg-promoter-methylation.Rmd')"

# Plot HGG H3 wt promoter methylation
R -e "rmarkdown::render('05-plot-gene-methylation.Rmd')"