# Assess germline PLP variant enrichment in OpenPBTA germline cohort

This module calculates enrichment of P-LP variants in cancer predisposition genes within PBTA germline cohort relative to two cancer-free control cohorts (Penn Med Biobank [PMBB] and gnomAD), and plots significant results. 

## Usage

`bash run_module.sh`

## Folder content 

1. `01-cpg-list-enr.R` plots CPG P-LP carrier enrichment in entire PBTA cohort and individual plot groups relative to tumor-free control cohorts

2. `02-gene-pathway-enrichment.R` plots gene- and pathway-level P-LP carrier enrichment in PBTA cohort relative to tumor-free control cohorts
### Input files

All files in `input/` directory were generated outside of this github repository, and utilize PBTA P-LP variant files to calculate enrichment of P-LP variants in the PBTA cohort relative to PMBB and gnomAD tumor-free control cohorts in the following gene lists:

1. All cancer predisposition genes (CPGs)
2. Individual CPGs ()
3. KEGG pathway gene sets
4. DNA repair pathway gene lists as reported in [Knijnenburg et al. 2018](https://www.cell.com/cell-reports/pdf/S2211-1247(18)30437-6.pdf)

3. `03-hist-gene-enr.R` plots gene-level P-LP carrier enrichment by histology cohort relative to tumor-free control cohorts

##Analysis module directory structure

```
.
├── 01-cpg-list-enr.R
├── 02-gene-pathway-enrichment.R
├── 03-hist-gene-enr.R
├── README.md
├── Rplots.pdf
├── input
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_cpg_pathway_pmbb_enrichment.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_dna_repair_pathway_pmbb_enrichment.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_gene_pmbb_enrichment.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-all-exome-filtered-20bp_padded_kegg_pathway_pmbb_enrichment.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-no-wxs_cpg_pathway_gnomad_enrichment.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-no-wxs_dna_repair_pathway_gnomad_enrichment.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged-no-wxs_gene_gnomad_enrichment.tsv
│   └── pbta-merged-plp-variants-autogvp-abridged-no-wxs_kegg_pathway_gnomad_enrichment.tsv
├── plots
│   ├── all-CPG-enrichment-PBTA-vs-control.tiff
│   ├── hist-all-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Atypical Teratoid Rhabdoid Tumor-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Choroid plexus tumor-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-KEGG_pathways-enrichment-PBTA-vs-control.tiff
│   ├── sig-Knijnenburg_repair_pathways-enrichment-PBTA-vs-control.tiff
│   ├── sig-Low-grade glioma-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Medulloblastoma-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Meningioma-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Mixed neuronal-glial tumor-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Neurofibroma plexiform-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Non-neoplastic tumor-CPG-enrichment-PBTA-vs-control.tiff
│   ├── sig-Other high-grade glioma-CPG-enrichment-PBTA-vs-control.tiff
│   └── sig-Pineoblastoma-CPG-enrichment-PBTA-vs-control.tiff
├── run_module.sh
└── util
    └── enrichment_functions.R
```