# Identify somatic second hits in PBTA germline cohort

This module queries SNV, CNV, LOH, and gene expression data from matched tumor samples of PBTA germline cohort patients to identify putative second hits in genes exhibiting germline pathogenic/likely pathogenic (P/LP) variation. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-assess-two-hits-snv.Rmd` identifies patients exhibiting both germline P/LP variation as well as an oncogenic/likely oncogenic somatic SNV in the same cancer predisposition gene (CPG). 

2. `02-cnv-loh-second-hits.Rmd` identifies patients exhibiting both germline P/LP variation as well as copy number variation (CNV) or loss of heterozygosity (LOH) in the same CPG. 

3. `input/` directory contains the following files: 
  - `PBTA-838-germline_plp-2023-02-05.tsv`; full germline P/LP results in PBTA cohort

4. `results/` directory contains the following files: 
  - `pbta-oncokb-oncogenic-maf.tsv`; subset of oncoKB results containing only oncogenic/likely oncogenic SNVs. 
  - `germline-somatic-two-gene-hits.tsv`; table of patients possessing germline P/LP variant and oncogenic/likely oncogenic SNV in same CPG. 
  - `germline-somatic-collapsed-by-gene.tsv`; summary of oncogenic/likely oncogenic SNVs in CPGs by patient and matched tumor sample. 
  - `germline-somatic-cnv-loh.tsv`; summary of CNV and LOH status by patients and CPG with germline P/LP variant. 
  
5. `plots/` directory contains the following files: 
  - `cpg-LOH-plp-vs-noplp.pdf`; box plots of LOH scores by P/LP status for select CPGs.
  - `hist-gene-LOH-plp-vs-noplp.pdf`; box plots of LOH scores by P/LP status for select CPGs within histologies.


##Analysis module directory structure

```
.
├── 01-assess-two-hits-snv.Rmd
├── 01-assess-two-hits-snv.nb.html
├── 02-cnv-loh-second-hits.Rmd
├── 02-cnv-loh-second-hits.nb.html
├── README.md
├── input
│   └── PBTA-838-germline_plp-2023-02-05.tsv
├── plots
│   ├── cpg-LOH-plp-vs-noplp.pdf
│   └── hist-gene-LOH-plp-vs-noplp.pdf
└── results
    ├── germline-somatic-cnv-loh-expr.tsv
    ├── germline-somatic-cnv-loh.tsv
    ├── germline-somatic-collapsed-by-gene.tsv
    ├── germline-somatic-two-gene-hits.tsv
    └── pbta-oncokb-oncogenic-maf.tsv
```
