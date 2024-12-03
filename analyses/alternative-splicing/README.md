# Identify germline P/LP variant-associated alternative splicing events in PBTA germline cohort

This module queries alternative splicing event RMATs file to identify germline P/LP variant-proximal alternative splicing events in matched tumors, and assesses associations between splice variants, alterantive splicing levels, and gene expression.

## Usage

`bash run_module.sh`

## Folder contents

1. `01-calculate-PSI-zscores.Rmd` calculates PSI z-scores of alternative splicing events proximal to germline P/LP variants in P/LP carriers

##Analysis module directory structure

```
.
├── 01-calculate-PSI-zscores.Rmd
├── 01-calculate-PSI-zscores.nb.html
├── README.md
├── input
│   └── splice_events_cpgs_pbta.tsv.gz
├── plots
│   ├── alternative-splicing-PSI-zscores-by-variant-type.pdf
├── results
│   ├── plp-variant-proximal-splicing-event-psi.tsv
└── run_module.sh
```
