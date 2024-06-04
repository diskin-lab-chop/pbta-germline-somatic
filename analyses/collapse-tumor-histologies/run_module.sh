#!/bin/bash

set -e
set -o pipefail

# Run collapse tumor histologies
Rscript 01-collapse-tumor-histologies.R

# Generate cohort summary plots
Rscript 02-cohort-summary-plots.R
