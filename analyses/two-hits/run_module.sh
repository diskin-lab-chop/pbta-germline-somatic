#!/bin/bash

set -e
set -o pipefail

# Run snv second hits analysis
Rscript -e "rmarkdown::render('01-assess-two-hits-snv.Rmd')"

# Run cnv/loh second hits analysis
Rscript -e "rmarkdown::render('02-cnv-loh-second-hits.Rmd')"

# Run somatic alteration enrichment analysis
Rscript -e "rmarkdown::render('03-somatic-alteration-enrichment.Rmd')"

# Run gene expression analysis
Rscript -e "rmarkdown::render('04-assess-somatic-gene-expr.Rmd')"