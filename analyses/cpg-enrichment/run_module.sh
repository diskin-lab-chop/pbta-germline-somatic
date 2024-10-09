#!/bin/bash

set -e
set -o pipefail

# Run CPG gene list enrichment
Rscript --vanilla 01-cpg-list-enr.R

# Run gene- and pathway-level enrichment
Rscript --vanilla 02-gene-pathway-enrichment.R