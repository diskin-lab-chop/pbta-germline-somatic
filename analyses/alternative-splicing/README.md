# Identify germline P/LP variant-associated alternative splicing events in PBTA germline cohort

This module queries alternative splicing event RMATs file to identify germline P/LP variant-proximal alternative splicing events in matched tumors, and assesses associations between splice variants, alterantive splicing levels, and gene expression.

## Usage

`bash run_module.sh`

## Folder contents

1. `01-calculate-PSI-zscores.Rmd` calculates PSI z-scores of alternative splicing events proximal to germline P/LP variants in P/LP carriers

2. `02-splicing-expression.Rmd` Plots gene expression TPM z-scores by variant classification and alternative splicing PSI z-scores, and identifies candidate P/LP variant-associated alternative splicing events. 

3. `03-plp-proximal-alternative-splicing.Rmd` identified significant P/LP variant-associated alternative splicing events, and calculates PSI-expression correlations in these events. 

##Analysis module directory structure

```
.
├── 01-calculate-PSI-zscores.Rmd
├── 01-calculate-PSI-zscores.nb.html
├── 02-splicing-expression.Rmd
├── 02-splicing-expression.nb.html
├── 03-plp-proximal-alternative-splicing.Rmd
├── 03-plp-proximal-alternative-splicing.nb.html
├── README.md
├── input
│   ├── gtf-annotated-splice-events.tsv
│   ├── splice_event_exon_intron_numbers.tsv
│   └── splice_events_cpgs_pbta.tsv.gz
├── plots
│   ├── ATM-exon 41-Low-grade glioma-tpm.pdf
│   ├── BAP1-exon 11-Meningioma-tpm.pdf
│   ├── FAH-exon 12-Low-grade glioma-tpm.pdf
│   ├── FAH-exon 12-Schwannoma-tpm.pdf
│   ├── FAH-exon 12-prot-abundance.pdf
│   ├── GBA-exon 2-Low-grade glioma-tpm.pdf
│   ├── GBA-exon 2-MB, SHH-tpm.pdf
│   ├── GBA-exon 2-prot-abundance.pdf
│   ├── LZTR1-exon 4-Atypical Teratoid Rhabdoid Tumor-tpm.pdf
│   ├── LZTR1-intron 14-Schwannoma-tpm.pdf
│   ├── NF1-exon 31-prot-abundance.pdf
│   ├── NF1-exon 5-Non-neoplastic tumor-tpm.pdf
│   ├── PTCH1-exon 2-MB, SHH-tpm.pdf
│   ├── SMARCE1-intron 8-9-Meningioma-tpm.pdf
│   ├── TSC2-exon 11-SEGA-tpm.pdf
│   ├── TSC2-exon 11-prot-abundance.pdf
│   ├── TSC2-exon 13-SEGA-tpm.pdf
│   ├── TSC2-exon 13-prot-abundance.pdf
│   ├── TSC2-exon 27-Non-neoplastic tumor-tpm.pdf
│   ├── TSC2-exon 32-SEGA-tpm.pdf
│   ├── TSC2-exon 32-prot-abundance.pdf
│   ├── TSC2-intron 11-SEGA-tpm.pdf
│   ├── TSC2-intron 11-prot-abundance.pdf
│   ├── TSC2-intron 26-27-Non-neoplastic tumor-tpm.pdf
│   ├── alternative-splicing-PSI-zscores-by-variant-type.pdf
│   ├── germline_plp_alternative_splicing_diff_plot.pdf
│   ├── tpm-altss-psi-correlations-by-variant-type.pdf
│   ├── tpm-psi-correlations-by-splicing-case.pdf
│   ├── tpm-retained-intron-psi-correlations-by-variant-type.pdf
│   ├── tpm-single-exon-psi-correlations-by-variant-type.pdf
│   └── tpm-zscores-by-variant-classification.pdf
├── results
│   ├── as-event-proximal-nonsplice-plp-variants.tsv
│   ├── candidate_plp_associated_splice_events.tsv
│   └── plp-variant-proximal-splicing-event-psi.tsv
├── run_module.sh
└── util
    └── splicing_plot_functions.R
```
