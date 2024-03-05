#!/bin/bash

set -e
set -o pipefail

# Summarize germline P/LP SVs
Rscript -e "rmarkdown::render('01-sv-summary.Rmd')"
