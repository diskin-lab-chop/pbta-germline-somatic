#!/bin/bash

set -e
set -o pipefail

# Run oncoprint script
Rscript -e "rmarkdown::render('01-oncoprint.Rmd')"
