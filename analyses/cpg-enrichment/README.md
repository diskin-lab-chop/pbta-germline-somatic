# Assess germline PLP variant enrichment in OpenPBTA germline cohort

This module calculates enrichment of P-LP variants in cancer predisposition genes within PBTA germline cohort relative to two cancer-free control cohorts (Penn Med Biobank and gnomAD), and plots significant results. 

## Usage

`bash run_module.sh`

## Folder content 

1. `01-cpg-list-enr.R` plots CPG P-LP carrier enrichment in entire PBTA cohort and individual plot groups relative to tumor-free control cohorts

2. `02-gene-pathway-enrichment.R` plots gene- and pathway-level P-LP carrier enrichment in PBTA cohort relative to tumor-free control cohorts


##Analysis module directory structure

```
.
├── 01-cpg-list-enr.R
├── 02-gene-pathway-enrichment.R
├── README.md
├── input
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_Knijnenburg_PMBB_enrichment.tsv
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_Knijnenburg_gnomAD_enrichment.tsv
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_PMBB_enrichment.tsv
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_gnomAD_enrichment.tsv
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_pathway_PMBB_enrichment.tsv
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_cpg_pathway_gnomAD_enrichment.tsv
│   ├── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_kegg_PMBB_enrichment.tsv
│   └── pbta-germline-837-plp-variants-nonpass-filtered-plus-reviewed_kegg_gnomAD_enrichment.tsv
├── plots
│   ├── all-CPG-enrichment-PBTA-vs-control.tiff
│   ├── hist-all-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-KEGG_pathways-enrichment-PBTA-vs-control.tiff
│   └── sig-Knijnenburg_repair_pathways-enrichment-PBTA-vs-control.tiff
├── run_module.sh
└── util
    └── enrichment_functions.R
```