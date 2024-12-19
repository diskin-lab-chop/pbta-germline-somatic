#!/bin/bash

set -e
set -o pipefail

# Format methylation data
Rscript --vanilla 01-prepare-methylation.R
