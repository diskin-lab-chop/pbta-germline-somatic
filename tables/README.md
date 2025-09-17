# Write manuscript tables

Module authors: Ryan Corbett (@rjcorb)

The purpose of this module is write and generate suppl tables for manuscript

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
Rscript --vanilla make-supp-tables.R
```

| Table | Supp file | Description | Location |
|-------|-----------|-------------|----------|
| TableS1_predisposition_syndromes     | NA        | CNS tumor predisposition syndromes defined by WHO 2021 Classification | NA (Manually curated) |
| TableS2_cohort-summary   | germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv        | Pediatric CNS tumor patient clinical and demographic data        | analyses/survival/results/ |
| TableS3_cpgs    | TableS3.xlsx        | Cancer predisposition gene chromosomal locations, associated syndromes and cancer types, and modes of inheritance  | tables/output/ |
| TableS4_plp-variant-summary    | pbta-merged-plp-variants-autogvp-abridged.tsv, pbta_germline_svs.tsv        | Germline P/LP variant summary  | data/v11/ |
| TableS5_syndrome_variants     | TableS5.xlsx        | Summary of predisposition syndrome-associated gene P/LP variants | tables/output/ |
| TableS6_syndrome_pts_VUS     | NA        | Germline variants of uncertain significance (VUS) and low variant allele frequency (VAF) P/LP variants in patients with clinically-reported syndromes but no associated P/LP variant | NA (manually curated) |
| TableS7_low_VAF_variants     | pbta-merged-plp-variants-autogvp-abridged-lowVAF-cpgs.tsv        | AutoGVP results of germline P/LP SNVs/indels of low allele frequency | data/v11/ |
| TableS8_plp_variants_by_molecular_subtype     | TableS8.xlsx        | CPG P/LP variants identified P/LP carrier-enriched molecular subtypes | tables/output/ |
| TableS9_demo-clin-stats-by-plp-carrier-status     | TableS9.xlsx        | Demographic and clinical features by CPG P/LP carrier status and histology | tables/output/ |
| TableS10_plp-variant-burden-testing     | TableS10.xlsx        | P/LP variant burden testing in PBTA versus Penn Medicine BioBank (PMBB) and gnomAD cancer-free control cohorts  | tables/output/ |
| TableS11_somatic_snv_indel_cnv_loh_expr    | TableS11.xlsx        | Summary of gene-level somatic alterations associated with germline P/LP variation | tables/output/ |
| TableS12_differential_methylation     | TableS12.xlsx        | Summary of germline CPG P/LP variant-associated differentially methylated Infinium 850K EPIC array probes in matched tumors  | tables/output/ |
| TableS13_alternative_splicing     | TableS13.xlsx        | Summary of germline CPG P/LP variant-proximal splice events in matched tumors | tables/output/ |
| TableS14_mutational_signatures     | TableS14.xlsx        | COSMIC single base substitution (SBS) mutational signatures exposures in HGG tumors of mismatch repair gene P/LP carriers versus non-carriers  | tables/output/ |
| TableS15_tmb_vs_methylation    | TableS15.xlsx        | Summary of tumor mutation – mean methylation beta value correlations by pediatric CNS tumor histology  | tables/output/ |
| TableS16_median-survival-by-plp-status     | TableS16.xlsx        | Median event- free survival (EFS) and overall survival (OS) in germline CPG P/LP carriers versus non-carriers by tumor histology and molecular subtype | tables/output/ |


## Folder content
* `make-suppl-tables.R` script to generate suppl tables in xls format for manuscript

## Directory structure
```
.
├── input
│   └── sbs_v3_map.tsv
├── make-supp-tables.R
├── output
│   ├── TableS10.xlsx
│   ├── TableS11.xlsx
│   ├── TableS12.xlsx
│   ├── TableS13.xlsx
│   ├── TableS14.xlsx
│   ├── TableS15.xlsx
│   ├── TableS16.xlsx
│   ├── TableS3.xlsx
│   ├── TableS5.xlsx
│   ├── TableS8.xlsx
│   └── TableS9.xlsx
├── README.html
└── README.md
```
