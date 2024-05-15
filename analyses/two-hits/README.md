# Identify somatic second hits in PBTA germline cohort

This module queries SNV, CNV, LOH, and gene expression data from matched tumor samples of PBTA germline cohort patients to identify putative second hits in genes exhibiting germline pathogenic/likely pathogenic (P/LP) variation. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-assess-two-hits-snv.Rmd` identifies patients exhibiting both germline P/LP variation as well as an oncogenic/likely oncogenic somatic SNV in the same cancer predisposition gene (CPG). 

2. `02-cnv-loh-second-hits.Rmd` identifies patients exhibiting both germline P/LP variation as well as copy number variation (CNV) or loss of heterozygosity (LOH) in the same CPG. 

4.`03-assess-somatic-gene-expr.Rmd`; pulls somatic gene expression data from matched tumor RNA-seq, and plots expression z-scores by P/LP status for select CPGs. 

5. `04-assess-alternative-splicing.Rmd`; pulls somatic splicing events and PSI values from matched tumor RNA-seq, and identifies alternative splicing events associated with germline CPG P/LP variants. 

6. `06-proteomics.Rmd`; pulls proteomics and phospoproteomics data from CPTAC and HOPE data sets, and identifies differential protein and phosphoprotein expression associated with germline CPG P/LP variants. 

7. `input/` directory contains the following files: 
  - `splice_events_germline_pbta_cpgs_only.tsv.gz`; splice events and PSI values for CPGs in PBTA germline samples
  - `HOPE_Harmonized_G_H_P_Proteome_UMich_RefSeq_Sinai_imputed_02-22-2022_V2.tsv`; HOPE proteomics total protein abundances
  - `PBT_Harmonized_G_H_P_Proteome_UMich_RefSeq_Sinai_imputed_06302022_exclude_bridging_V2.tsv`; CPTAC proteomics total protein abundances
  - `hope-protein-imputed-phospho-expression.tsv`; HOPE total phosphoprotein abundances 
  - `histologies-base.tsv`; v13 histologies base file (will be removed here following v7 data release)

8. `results/` directory contains the following files: 
  - `pbta-oncokb-oncogenic-maf.tsv`; subset of oncoKB results containing only oncogenic/likely oncogenic SNVs
  - `germline-somatic-two-gene-hits.tsv`; table of patients possessing germline P/LP variant and oncogenic/likely oncogenic SNV in same CPG
  - `germline-somatic-collapsed-by-gene.tsv`; summary of oncogenic/likely oncogenic SNVs in CPGs by patient and matched tumor sample
  - `germline-somatic-cnv-loh.tsv`; summary of CNV and LOH status by patients and CPG with germline P/LP variant
  - `somatic-alteration-enrichment.tsv`; summary of somatic alteration enrichment scores (Fisher's exact tests) by gene and alteration type
  - `pbta-germline-838-gene-expr-zscores.tsv`; expression z-scores in PBTA germline cohort
  - `germline-somatic-expr.tsv`; summary of expression z-scores for patient and CPG harboring a germline P/LP variant
  - `gene-expr-diff-plp-vs-no-plp.tsv`; mean expression z-scores by gene and germline P/LP status with wilcoxon rank p-values assessing significance of difference  
  - `splicing_events_plp_variants.tsv`; splicing events and PSI values, z-scores, and differences relative to other hist group samples in in CPG P/LP carriers
  - `candidate_plp_associated_splice_events.tsv`; plp variants proximal (<250bp) from differential splice event
  - `cptac-hope-cpg-proteomics-zscores.tsv`; protein abundance z-scores in PBTA germline cohort
  - `hope-cpg-phosphoproteomics-zscores.tsv`; phosphoprotein abunddance z-scores in PBTA germline cohort
  - `germline-somatic-proteomics-phosphoproteomics-cpgs.tsv`; protein and phosphoprotein abundance z-scores in patietns with germline P/LP variant in same CPG. 
  
  
9. `plots/` directory contains the following files: 
  - `cpg-LOH-plp-vs-noplp.pdf`; box plots of LOH scores by P/LP status for select CPGs.
  - `hist-gene-LOH-plp-vs-noplp.pdf`; box plots of LOH scores by P/LP status for select CPGs within histologies.
  - `cpg-LOH-plp-vs-noplp.pdf`; box plots of LOH score by germline P/LP status for select CPGs.
  - `hist-gene-LOH-plp-vs-noplp.pdf` box plots of LOH score by germline P/LP status for select CPGs within histologies. 
  - `cpg-sig-expr-diff-plp-vs-noplp.pdf`; box plots of expression z-scores by germline P/LP status for select CPGs. 
  - `hist-gene-expr-plp-vs-noplp.pdf`; box plots of expression z-scores by germline P/LP status for select CPGs within histologies.
  - `*-expr-by-plp-onco-snv-status.pdf`; gene-specific box plots of expression z-scores by germline P/LP status and somatic oncogenic/likely oncogenic SNV status. 
  - `PSI-diff-by-variant-dist.pdf`; scatter plot of PSI difference by P/LP variant distance to splice region. 


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
├── README.md
├── input
│   ├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
│   ├── histologies-base.tsv
│   ├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
│   ├── hope-protein-imputed-prot-expression-abundance.tsv.gz
│   └── splice_events_germline_pbta_cpgs_only.tsv.gz
├── plots
│   ├── ATM-expr-by-plp-onco-snv-status.pdf
│   ├── CHEK2-expr-by-plp-onco-snv-status.pdf
│   ├── ERCC2-expr-by-plp-onco-snv-status.pdf
│   ├── FANCA-expr-by-plp-onco-snv-status.pdf
│   ├── MUTYH-expr-by-plp-onco-snv-status.pdf
│   ├── NF1-expr-by-plp-onco-snv-status.pdf
│   ├── NF2-expr-by-plp-onco-snv-status.pdf
│   ├── PMS2-expr-by-plp-onco-snv-status.pdf
│   ├── PSI-diff-by-variant-dist.pdf
│   ├── PSI-diff-by-variant-dist_all_plp.pdf
│   ├── RAD50-expr-by-plp-onco-snv-status.pdf
│   ├── RECQL4-expr-by-plp-onco-snv-status.pdf
│   ├── TP53-expr-by-plp-onco-snv-status.pdf
│   ├── TSC1-expr-by-plp-onco-snv-status.pdf
│   ├── TSC2-expr-by-plp-onco-snv-status.pdf
│   ├── cpg-LOH-plp-vs-noplp.pdf
│   ├── cpg-sig-expr-diff-plp-vs-noplp.pdf
│   ├── hist-gene-LOH-plp-vs-noplp.pdf
│   ├── hist-gene-expr-plp-vs-noplp.pdf
│   ├── proximal-splice-event-distribution.tiff
│   ├── sig-somatic-alteration-enr.pdf
│   └── somatic-alteration-enr-heatmap.pdf
├── results
│   ├── candidate_plp_associated_splice_events.tsv
│   ├── cptac-hope-cpg-proteomics-zscores.tsv
│   ├── gene-expr-diff-plp-vs-no-plp.tsv
│   ├── germline-somatic-cnv-loh-cpgs.tsv
│   ├── germline-somatic-cnv-loh-somatic-non-cpgs.tsv
│   ├── germline-somatic-collapsed-by-gene.tsv
│   ├── germline-somatic-expr-cpgs.tsv
│   ├── germline-somatic-expr-somatic-non-cpgs.tsv
│   ├── germline-somatic-proteomics-phosphoproteomics-cpgs.tsv
│   ├── germline-somatic-two-gene-hits.tsv
│   ├── hope-cpg-phosphoproteomics-zscores.tsv
│   ├── pbta-germline-838-gene-expr-zscores.tsv
│   ├── pbta-oncokb-oncogenic-maf.tsv
│   └── splicing_events_plp_variants.tsv
└── run_module.sh
```
