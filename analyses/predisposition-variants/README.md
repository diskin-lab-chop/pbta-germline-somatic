# Assess cancer predisposition syndrome-associated gene P/LP variant frequency in PBTA germline cohort

This module identifies all P/LP variants in genes associated with cancer predisposition syndromes, and cross-references findings with clinical reports of syndromes to determine if findings are incidental. Conversely, the frequency of patients with a reported syndrome AND that have an identified P/LP variant in the associated gene is plotted. 

## Usage

`bash run_module` 

## Folder content 

1. `01-incidental-findings.Rmd`; pulls all variants in syndrome-associated genes, and determines if finding is incidental by integrating histology data

2. `02-plot-predispositions.Rmd`; creates predisposition variant plots

## Directory Structure
```
.
├── 01-incidental-findings.Rmd
├── 01-incidental-findings.html
├── 02-plot-predispositions.Rmd
├── 02-plot-predispositions.html
├── README.md
├── input
│   ├── predisposition-syndromes.tsv
│   └── predispositions-path-review.tsv
├── plots
│   ├── gene-plp-prevalence-by-syndrome.pdf
│   └── syndrome-prevalence-by-gene.pdf
├── results
│   ├── incidental-findings-predisposition-structural-variants-path-review.tsv
│   ├── incidental-findings-predisposition-structural-variants.tsv
│   ├── incidental-findings-predisposition-variants-path-review.tsv
│   ├── incidental-findings-predisposition-variants.tsv
│   └── predisposition-patients-no-plp-variants.tsv
└── run_module.sh
```