# PBTA Germline: methylation data integration

This module processes DNA methylation array data and tests for differential methylation between germline P-LP carriers and non-carriers by histology

## Usage

`bash run_module.sh` 

## Folder content 

1. `01-prepare_methylation.R` Calculate global and mean promoter and gene body methylation rates in CPGs, and extract all CPG beta values.

2. `02-global-methylation.Rmd` Compare mean beta values of P-LP carriers to non-carriers within histology groups

3. `03-calculate-probe-zscores.Rmd` Compare probe beta values of P-LP carriers to non-carriers within histology groups

## Directory structure

```
.
├── 01-prepare-methylation.R
├── 02-global-methylation.Rmd
├── 02-global-methylation.nb.html
├── 03-calculate-probe-zscores.Rmd
├── 03-calculate-probe-zscores.nb.html
├── plots
│   ├── global-beta-value-BER-plp-vs-other.pdf
│   ├── global-beta-value-DNA repair-plp-vs-other.pdf
│   ├── global-beta-value-HR-plp-vs-other.pdf
│   ├── global-beta-value-MMR-plp-vs-other.pdf
│   ├── global-beta-value-NER-plp-vs-other.pdf
│   ├── global-beta-value-NHEJ-plp-vs-other.pdf
│   ├── global-beta-value-hgg-h3wt-mmr-plp-vs-other.pdf
│   ├── mean-sample-methyl-beta-by-hist-gene-plp.pdf
│   ├── mean-sample-methyl-beta-by-histology-plp.pdf
│   ├── mean-sample-methyl-beta-by-histology.pdf
│   └── tmb-vs-mean-beta-methylation-hgg-h3-wt.pdf
├── results
│   ├── annot-with-canonical.tsv
│   ├── cpg-methyl-beta-values.rds
│   ├── gene-methyl-zscores.rds
│   ├── pbta-germline-mean-sample-methyl-beta.tsv
│   ├── pbta-germline-methylation-beta-zscores.tsv
│   └── promoter-methyl-zscores.rds
└── run_module.sh
```
