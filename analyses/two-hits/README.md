# Identify somatic second hits in PBTA germline cohort

This module queries SNV, CNV, LOH, and gene expression data from matched tumor samples of PBTA germline cohort patients to identify putative second hits in genes exhibiting germline pathogenic/likely pathogenic (P/LP) variation. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-assess-two-hits-snv.Rmd` identifies patients exhibiting both germline P/LP variation as well as an oncogenic/likely oncogenic somatic SNV in the same cancer predisposition gene (CPG). 

2. `02-cnv-loh-second-hits.Rmd` identifies patients exhibiting both germline P/LP variation as well as copy number variation (CNV) or loss of heterozygosity (LOH) in the same CPG. 

3. `03-somatic-alteration-enrichment.Rmd` tests for enrichment of somatic alterations (SNV, CNV, LOH, gene expression gain/loss) in patients with vs. without germline P/LP variants in select CPGs. 

4.`04-assess-somatic-gene-expr.Rmd`; pulls somatic gene expression data from matched tumor RNA-seq, and plots expression z-scores by P/LP status for select CPGs. 

5. `05-assess-alternative-splicing.Rmd`; pulls somatic splicing events and PSI values from matched tumor RNA-seq, and identifies alternative splicing events associated with germline CPG P/LP variants. 

6. `06-proteomics.Rmd`; pulls proteomics and phospoproteomics data from CPTAC and HOPE data sets, and identifies differential protein and phosphoprotein expression associated with germline CPG P/LP variants.

7. `07-merge-somatic-alterations.Rmd`; merges all P/LP variant and same-gene somatic alteration data for P/LP carriers

##Analysis module directory structure

```
.
├── 01-assess-two-hits-snv.Rmd
├── 01-assess-two-hits-snv.nb.html
├── 02-cnv-loh-second-hits.Rmd
├── 02-cnv-loh-second-hits.nb.html
├── 03-somatic-alteration-enrichment.Rmd
├── 03-somatic-alteration-enrichment.nb.html
├── 04-assess-somatic-gene-expr.Rmd
├── 04-assess-somatic-gene-expr.nb.html
├── 05-assess-alternative-splicing.Rmd
├── 05-assess-alternative-splicing.nb.html
├── 06-proteomics.Rmd
├── 06-proteomics.nb.html
├── 07-merge-somatic-alterations.Rmd
├── 07-merge-somatic-alterations.nb.html
├── README.md
├── input
│   ├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
│   ├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
│   ├── hope-protein-imputed-prot-expression-abundance.tsv.gz
│   ├── pbta-plp-all-loh.tsv
│   └── splice_events_germline_pbta_cpgs_only.tsv.gz
├── plots
│   ├── APC-expr-by-plp-onco-snv-status.pdf
│   ├── ATM-expr-by-plp-onco-snv-status.pdf
│   ├── CHEK2-expr-by-plp-onco-snv-status.pdf
│   ├── ERCC2-expr-by-plp-onco-snv-status.pdf
│   ├── FANCA-expr-by-plp-onco-snv-status.pdf
│   ├── MUTYH-expr-by-plp-onco-snv-status.pdf
│   ├── NF1-expr-by-plp-onco-snv-status.pdf
│   ├── NF2-expr-by-plp-onco-snv-status.pdf
│   ├── PMS2-expr-by-plp-onco-snv-status.pdf
│   ├── PSI-diff-by-variant-dist.pdf
│   ├── RAD50-expr-by-plp-onco-snv-status.pdf
│   ├── RECQL4-expr-by-plp-onco-snv-status.pdf
│   ├── TP53-expr-by-plp-onco-snv-status.pdf
│   ├── TSC1-expr-by-plp-onco-snv-status.pdf
│   ├── TSC2-expr-by-plp-onco-snv-status.pdf
│   ├── cpg-LOH-plp-vs-noplp.pdf
│   ├── cpg-sig-expr-diff-plp-vs-noplp.pdf
│   ├── hist-gene-LOH-plp-vs-noplp.pdf
│   ├── hist-gene-expr-plp-vs-noplp.pdf
│   ├── sig-somatic-alteration-enr.pdf
│   ├── skipped-exon-psi-splice-plp-vs-nonsplice.pdf
│   └── somatic-alteration-enr-heatmap.pdf
├── results
│   ├── candidate_plp_associated_splice_events.tsv
│   ├── candidate_sv_plp_associated_splice_events.tsv
│   ├── cptac-hope-cpg-proteomics-zscores.tsv
│   ├── gene-expr-diff-plp-vs-no-plp.tsv
│   ├── germline-primary-plus-tumor-histologies-proteomics.tsv
│   ├── germline-somatic-cnv-loh-cpgs.tsv
│   ├── germline-somatic-cnv-loh-somatic-non-cpgs.tsv
│   ├── germline-somatic-collapsed-by-gene.tsv
│   ├── germline-somatic-expr-cpgs.tsv
│   ├── germline-somatic-expr-somatic-non-cpgs.tsv
│   ├── germline-somatic-proteomics-cpgs.tsv
│   ├── germline-somatic-proteomics-phosphoproteomics-cpgs.tsv
│   ├── germline-somatic-two-gene-hits.tsv
│   ├── germline-sv-somatic-cnv-loh-somatic-cpgs.tsv
│   ├── germline-sv-somatic-expr-cpgs.tsv
│   ├── germline-sv-somatic-proteomics-phosphoproteomics-cpgs.tsv
│   ├── germline-sv-somatic-two-gene-hits.tsv
│   ├── hope-cpg-phosphoproteomics-zscores.tsv
│   ├── pbta-germline-838-gene-expr-zscores.tsv
│   ├── pbta-oncokb-oncogenic-maf.tsv
│   ├── plp_snv_indel_somatic_alterations_merged.tsv
│   ├── plp_sv_somatic_alterations_merged.tsv
│   ├── somatic-alteration-enrichment.tsv
│   ├── splicing_events_plp_variants.tsv
│   ├── splicing_events_snv_plp_variants.tsv
│   └── splicing_events_sv_plp_variants.tsv
└── run_module.sh
```
