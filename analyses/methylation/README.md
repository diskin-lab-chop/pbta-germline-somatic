# PBTA Germline: methylation data integration

This module processes DNA methylation array data and tests for differential methylation between germline P-LP carriers and non-carriers by histology

## Usage

`bash run_module.sh` 

## Folder content 

1. `01-prepare_methylation.Rmd` Calculate global and mean promoter and gene body methylation rates in CPGs, and extract all CPG beta values.

## Directory structure

```
.
├── 01-prepare-methylation.R
├── plots
├── results
│   ├── annot-with-canonical.tsv
│   ├── cpg-methyl-beta-values.rds
│   ├── gene-methyl-zscores.rds
│   ├── pbta-germline-mean-sample-methyl-beta.tsv
│   └── promoter-methyl-zscores.rds
└── run_module.sh
```
