#!/bin/bash

set -e
set -o pipefail

# set input directory
input_dir=input

awk '{start=$2-20; if (start < 0) start=0; print "chr"$1"\t"start"\t"$3+20}' $input_dir/xgen_plus_spikein.b38.bed > $input_dir/xgen_plus_spikein.b38.20bp_padded.bed
awk '{start=$2-100; if (start < 0) start=0; print "chr"$1"\t"start"\t"$3+100}' $input_dir/xgen_plus_spikein.b38.bed > $input_dir/xgen_plus_spikein.b38.100bp_padded.bed
