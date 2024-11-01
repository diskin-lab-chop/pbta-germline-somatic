#!/bin/bash

set -e
set -o pipefail

# Run CPG gene list enrichment
Rscript --vanilla 01-cpg-list-enr.R

# Run gene- and pathway-level enrichment
Rscript --vanilla 02-gene-pathway-enrichment.R

# Run gene-level enrichment by histology
Rscript --vanilla 03-hist-gene-enr.R

# Run pathway-level enrichment by histology
Rscript --vanilla 04-hist-pathway-enr.R