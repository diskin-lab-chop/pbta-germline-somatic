#!/bin/bash

set -e
set -o pipefail

# Format survival data
Rscript --vanilla 01-prepare_survival.R
