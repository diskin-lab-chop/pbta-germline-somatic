#!/bin/bash

set -e
set -o pipefail

# Get exome capture region intersection
bash 01-run-pad-intersect-beds.sh

# filter WXS P-LP variants
R -e "rmarkdown::render('02-filter-wxs-variants.Rmd')"