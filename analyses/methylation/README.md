# PBTA Germline: methylation data integration

This module processes DNA methylation array data and tests for differential methylation between germline P-LP carriers and non-carriers by histology

## Usage

`bash run_module.sh` 

## Folder content 

1. `00-rename-methylation-columns.R` One-time preprocessing step to rename columns within the methylation data file for v11 data release.
2. `01-prepare_methylation.R` Calculate global and mean promoter and gene body methylation rates in CPGs, and extract all CPG beta values.
3. `02-global-methylation.Rmd` Compare mean beta values of P-LP carriers to non-carriers within histology groups
4. `03-calculate-probe-zscores.Rmd` Compare probe beta values of P-LP carriers to non-carriers within histology groups
5. `04-plot-hgg-promoter-methylation.Rmd` Generates a heatmap of HGG H3 Wt samples by promoter methylation
6. `05-plot-gene-methylation.Rmd`; Plot regions of P/LP variant-associated differential methylation

## Directory structure

```
.
├── 00-rename-methylation-columns.R
├── 01-prepare-methylation.R
├── 02-global-methylation.nb.html
├── 02-global-methylation.Rmd
├── 03-calculate-probe-zscores.nb.html
├── 03-calculate-probe-zscores.Rmd
├── 04-plot-hgg-promoter-methylation.nb.html
├── 04-plot-hgg-promoter-methylation.Rmd
├── 05-plot-gene-methylation.nb.html
├── 05-plot-gene-methylation.Rmd
├── input
│   └── itt-141-pbta-methylation.tsv
├── plot-gene-methylation.nb.html
├── plot-gene-methylation.Rmd
├── plots
│   ├── expr-zscore-by-diff-meth.pdf
│   ├── global-beta-value-BER-plp-vs-other.pdf
│   ├── global-beta-value-DNA repair-plp-vs-other.pdf
│   ├── global-beta-value-hgg-h3wt-mmr-plp-vs-other.pdf
│   ├── global-beta-value-HR-plp-vs-other.pdf
│   ├── global-beta-value-MMR-plp-vs-other.pdf
│   ├── global-beta-value-NER-plp-vs-other.pdf
│   ├── global-beta-value-NHEJ-plp-vs-other.pdf
│   ├── mean-sample-methyl-beta-by-hist-gene-plp.pdf
│   ├── mean-sample-methyl-beta-by-histology-plp.pdf
│   ├── mean-sample-methyl-beta-by-histology.pdf
│   ├── MSH2-BS_RJMW12ZA-promoter-methylation.pdf
│   ├── msh2-gene-model.pdf
│   ├── MSH6-BS_5V730EQB-promoter-methylation.pdf
│   ├── msh6-gene-model-BS_G113MX63.pdf
│   ├── pbta-germline-hgg-wt-plp-somatic-prom-methylation.pdf
│   ├── plp-diff-probes-by-feature.pdf
│   ├── RECQL4-BS_KGT59ZDT-promoter-methylation.pdf
│   ├── RECQL4-BS_RWKHJFYY-promoter-methylation.pdf
│   ├── recql4-gene-model-BS_PF7546KW.pdf
│   ├── recql4-gene-model-BS_TV8HD883.pdf
│   ├── tmb-vs-mean-beta-methylation-hgg-h3-mutant.pdf
│   └── tmb-vs-mean-beta-methylation-hgg-h3-wt.pdf
├── README.md
├── results
│   ├── annot-with-canonical.tsv
│   ├── cpg-methyl-beta-values.rds
│   ├── gene-methyl-zscores.rds
│   ├── germline-primary-plus-tumor-histologies-methylation.tsv
│   ├── pbta-germline-mean-sample-methyl-beta.tsv
│   ├── pbta-germline-methylation-beta-zscores-plus-expr.tsv
│   ├── pbta-germline-methylation-beta-zscores.tsv
│   ├── promoter-methyl-zscores.rds
│   └── tmb-methylation-correlation-by-hist.tsv
└── run_module.sh
```
