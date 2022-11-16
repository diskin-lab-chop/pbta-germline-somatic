#!/bin/bash

set -e
set -o pipefail

python_path=(/usr/bin/python3)
maf_annotator=(/home/oncokb-annotator/MafAnnotator.py)
maf_in=(results/snv-consensus-plus-hotspots-goi.maf.tsv)
maf_oncokb_out=(results/snv-consensus-plus-hotspots-goi-oncokb.maf.tsv)

#OncoKB access token required; if not provided then exit 
if [ "$ONCOKB" = "" ]; then
    echo 'Error: An Oncokb access token is required to run this script. Please provide token before running script in the format: ONCOKB=<TOKEN> bash run-oncokb-annotation.sh'
fi

# Run maf_annotator on dgd samples
$python_path $maf_annotator -i $maf_in -o $maf_oncokb_out -b $ONCOKB -q hgvsp_short -r GRCh38
