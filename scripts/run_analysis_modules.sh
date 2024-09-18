#!/bin/sh

set -e
set -o pipefail

# OncoKB access token required to run module; if not provided then exit
if [ "$ONCOKB" = "" ]; then
    echo 'Error: An Oncokb access token is required to run this module in full. Please provide token before running analysis module in the format ONCOKB=<TOKEN> bash run-oncokb-annotation.sh'
    exit 1
fi

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

# Run collapse-tumor-histologies analysis module
echo "Run collapse-tumor-histologies"
cd ${analyses_dir}/collapse-tumor-histologies
Rscript 01-collapse-tumor-histologies.R

# Run bed-intersect analysis module
echo "Run bed-intersect"
cd ${analyses_dir}/bed-intersect
bash run_module.sh

# Run variant-distribution analysis module
echo "Run variant-distribution"
cd ${analyses_dir}/variant-distribution
bash run_module.sh

# Run germline-sv analysis module
echo "Run germline-sv"
cd ${analyses_dir}/germline-sv
bash run_module.sh

# Run sample-distribution analysis module
echo "Run sample-distribution"
cd ${analyses_dir}/sample-distribution
bash run_module.sh

# Run predisposition variant analysis module
echo "Run predisposition-variants"
cd ${analyses_dir}/predisposition-variants
bash run_module.sh

# Run enrichment analysis module
# echo "Run cpg-enrichment"
# cd ${analyses_dir}/cpg-enrichment
# bash run_module.sh

# Run oncokb annotation analysis module
echo "Run oncokb-annotation"
cd ${analyses_dir}/oncokb-annotation
ONCOKB=$ONCOKB bash run-oncokb-annotation.sh

# Run two-hits analysis module
echo "Run two-hits"
cd ${analyses_dir}/two-hits
bash run_module.sh

# Run two-hits analysis module
echo "Run oncoprint"
cd ${analyses_dir}/oncoprint
bash run_module.sh

# Run DNA repair variant variant analysis module
echo "Run dna-repair-variant-summary"
cd ${analyses_dir}/dna-repair-variant-summary
bash run_dna_repair_variant_analysis.sh

# Run survival analysis module
echo "Run survival"
cd ${analyses_dir}/survival
bash run_module.sh

# Run demo-clin stats analysis module
echo "Run demo clin stats"
cd ${analyses_dir}/demo-clin-stats
bash run_module.sh

printf "\nDone running analysis modules!\n"