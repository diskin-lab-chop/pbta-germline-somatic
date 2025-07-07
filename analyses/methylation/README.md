# PBTA Germline: methylation data integration

This module processes DNA methylation array data and tests for differential methylation between germline P-LP carriers and non-carriers by histology

## Usage

`bash run_module.sh` 

## Folder content 

1. `00-rename-methylation-columns.R` One-time preprocessing step to rename columns within the methylation data file for v11 data release.

2. `01-prepare_methylation.R` Calculate global and mean promoter and gene body methylation rates in CPGs, and extract all CPG beta values.

3. `02-global-methylation.Rmd` Compare mean beta values of P-LP carriers to non-carriers within histology groups

4. `03-calculate-probe-zscores.Rmd` Compare probe beta values of P-LP carriers to non-carriers within histology groups

4. `04-plot-hgg-promoter-methylation.Rmd` Generates a heatmap of HGG H3 Wt samples by promoter methylation

## Directory structure

```
.
├── 01-prepare-methylation.R
├── 02-global-methylation.Rmd
├── 02-global-methylation.nb.html
├── 03-calculate-probe-zscores.Rmd
├── 03-calculate-probe-zscores.nb.html
├── 04-plot-hgg-promoter-methylation.Rmd
├── 04-plot-hgg-promoter-methylation.nb.html
├── README.md
├── input
│   └── itt-141-pbta-methylation.tsv
├── plots
│   ├── expr-zscore-by-diff-meth.pdf
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
│   ├── pbta-germline-hgg-wt-plp-somatic-prom-methylation.pdf
│   ├── plp-diff-probes-by-feature.pdf
│   ├── tmb-vs-mean-beta-methylation-hgg-h3-mutant.pdf
│   └── tmb-vs-mean-beta-methylation-hgg-h3-wt.pdf
├── results
│   ├── annot-with-canonical.tsv
│   ├── cpg-methyl-beta-values.rds
│   ├── gene-methyl-zscores.rds
│   ├── germline-primary-plus-tumor-histologies-methylation.tsv
│   ├── hgg-methyl-beta-value.tsv.gz
│   ├── pbta-germline-mean-sample-methyl-beta.tsv
│   ├── pbta-germline-methylation-beta-zscores-plus-expr.tsv
│   ├── pbta-germline-methylation-beta-zscores.tsv
│   ├── promoter-methyl-zscores.rds
│   └── tmb-methylation-correlation-by-hist.tsv
└── run_module.sh
```
