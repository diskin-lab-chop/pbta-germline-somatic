# Perform bedtools intersect on WXS region files

This module pads the PMBB BED calling regions and does an intersection of all bed files to get consensus calling regions for WXS samples.

## Usage

`bash pad-intersect-beds.sh` 

## Folder content 

1. `pad-intersect-beds.sh` pads the pmbb bed file with either 20 or 100 bp and performs an intersection of bed regions across the 4 bed files:


## Directory Structure
```
.
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
├── pad-intersect-beds.sh
└── results
```
  
