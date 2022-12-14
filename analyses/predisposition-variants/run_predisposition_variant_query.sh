#!/bin/bash

# J. Taroni for ALSF CCDL 
# 2020

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

#Query variants that pass filtering in cancer predisposition patients 
Rscript -e "rmarkdown::render('01-query-predisposition-variants.Rmd', clean = TRUE)"

#Query all unfiltered variants in cancer predisposition patients 
Rscript -e "rmarkdown::render('02-query-unfiltered-predisposition-variants.Rmd', clean = TRUE)"