#!/bin/bash

set -e
set -o pipefail

# set input directory
input_dir=input
results_dir=results

awk '{start=$2-20; if (start < 0) start=0; print "chr"$1"\t"start"\t"$3+20}' $input_dir/xgen_plus_spikein.b38.bed > $input_dir/xgen_plus_spikein.b38.20bp_padded.bed
awk '{start=$2-100; if (start < 0) start=0; print "chr"$1"\t"start"\t"$3+100}' $input_dir/xgen_plus_spikein.b38.bed > $input_dir/xgen_plus_spikein.b38.100bp_padded.bed

# intersect regions from all 4 bed files
bedtools intersect -a $input_dir/xgen_plus_spikein.b38.20bp_padded.bed -b $input_dir/xgen-exome_hg38_liftover.merged.100bp_padded.bed \
  | bedtools intersect -a - -b $input_dir/ashion_confidential_exome_v2_targets_hg38_paded100.bed \
  | bedtools intersect -a - -b $input_dir/S07604715_100bp_Padded.bed > $results_dir/pmbb_pbta_wxs_intersect.bed
