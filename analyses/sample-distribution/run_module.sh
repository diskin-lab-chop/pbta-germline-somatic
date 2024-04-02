#!/bin/bash

set -e
set -o pipefail

# Run sample distribution 
Rscript --vanilla 01-create-cancer-distribution-plot.R
