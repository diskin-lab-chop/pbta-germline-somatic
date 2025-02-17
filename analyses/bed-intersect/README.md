# Perform bedtools intersect on WXS region files and filters variants for capture region overlap 

This module pads the PMBB BED calling regions and does an intersection of all bed files to get consensus calling regions for WXS samples. P-LP variants are subsequently filtered for overlap with these consensus calling regions. 

## Usage

`bash run_module.sh` 

## Folder content 

1. `01-pad-intersect-beds.sh` pads the pmbb bed file with either 20 or 100 bp and performs an intersection of bed regions across the 4 bed files:

2. `02-filter-variants.Rmd` Performs intersection of P-LP variants and consensus calling regions and filters out any non-overlapping variants from P-LP variant file


## Directory Structure
```
.
├── 01-run-pad-intersect-beds.sh
├── 02-filter-variants.Rmd
├── 02-filter-variants.html
├── README.md
├── input
│   ├── S07604715_100bp_Padded.bed
│   ├── S07604715_Regions.bed
│   ├── ashion_confidential_exome_v2_targets_hg38.bed
│   ├── ashion_confidential_exome_v2_targets_hg38_paded100.bed
│   ├── region_40samples.csv
│   ├── xgen-exome-research-panel-targets_hg38_ucsc_liftover.sort.merged.bed
│   ├── xgen-exome_hg38_liftover.merged.100bp_padded.bed
│   ├── xgen_plus_spikein.b38.100bp_padded.bed
│   ├── xgen_plus_spikein.b38.20bp_padded.bed
│   └── xgen_plus_spikein.b38.bed
├── results
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-excluded-100bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-excluded-20bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-100bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-no-wxs.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-pbta-exome-excluded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-pbta-exome-filtered.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-wxs-exome-excluded-100bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-wxs-exome-excluded-20bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-wxs-exome-filtered-100bp_padded.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-wxs-exome-filtered-20bp_padded.tsv
│   ├── pbta_wxs_intersect.bed
│   ├── pmbb_pbta_wxs_intersect.100bp_padded.bed
│   └── pmbb_pbta_wxs_intersect.20bp_padded.bed
└── run_module.sh
```
  
