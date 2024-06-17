#!/bin/bash

set -e
set -o pipefail

# Run snv second hits analysis
Rscript -e "rmarkdown::render('01-assess-two-hits-snv.Rmd')"

# Run cnv/loh second hits analysis
Rscript -e "rmarkdown::render('02-cnv-loh-second-hits.Rmd')"

# Run gene expression analysis
Rscript -e "rmarkdown::render('03-assess-somatic-gene-expr.Rmd')"

# Run gene expression analysis
Rscript -e "rmarkdown::render('04-assess-alternative-splicing.Rmd')"

# Run proteomics analysis
Rscript -e "rmarkdown::render('05-proteomics.Rmd')"

# merge alterations output
Rscript -e "rmarkdown::render('06-merge-somatic-alterations.Rmd')"
