# Characterization of germline P-LP variation in DNA repair genes 

This module summarizes germline pathogenic variation in known DNA repair genes, and further compares mutational signatures and pathway expression in HGG tumors with vs. without P-LP variants in DNA repair pathway genes. 

## Usage

`bash run_dna_repair_variant_analysis.sh` 

## Folder content 

1. `01-get-dna-repair-variant-samples.Rmd` identifies and plots summary of germline samples with P-LP variants in DNA repair genes. 
2. `02-mutational-signatures.Rmd` compares mutational signature exposure levels in high-grade glioma patients with versus without germline P-LP variants in DNA repair genes.
3. `03-run-gsva-comparison.Rmd` compares gene set variation analysis (GSVA) scores in high-grade glioma patients with versus without germline P-LP variants in DNA repair genes.
4. `run-kegg-gsva.R` runs GSVA analysis of KEGG pathways in PBTA samples. 

4. `input/` directory contains the following files:
  - `COSMICv3.3_signature_exposures_SBS39_included.tsv`
  - `broad/` and `Knijnenburg_paper/` contain the following DNA repair pathway gene lists from GSEA and Knijnenburg et al. 2021:
      - `dna_repair_all.txt`
      - `base_excision_repair.txt`
      - `homologous_recombination.txt`
      - `mismatch_repair.txt`
      - `nonhomologous_end_joining.txt`
      - `nucleotide_excision_repair.txt`
  - `pbta-hallmark-gsva-scores.tsv`
    

3. `results/` directory contains the following files: in both `Broad_GO_repair_genes/` and `Knijnenberg_repair_genes/`: 
  - `pbta-kegg-gsva-scores.tsv`
  - The following files are found in both `Broad_GO_repair_genes/` and `Knijnenberg_repair_genes/` subdirectories:
    - mutational signature difference tables in HGG tumors of patients with vs without DNA repair gene germline P-LP variants. 
    - gsva score difference tables in HGG tumors of patients with vs without DNA repair gene germline P-LP variants. 
  
4. `plots/` directory contains the following figures in both `Broad_GO_repair_genes/` and `Knijnenberg_repair_genes/`: 
  - Summary and enrichment of DNA repair genes by histology (`dna_repair_samples_by_plot_group_*` and `dna_repair_enr_by_plot_group_*`)
  - Bar plots of COSMICv3 signature exposure differences in HGG patients with versus without P-LP variants in DNA repair genes. 
  - Violin plots of COSMICv3 signature exposure weights in HGG patients with 1) MMR gene P-LP variants, 2) BRCA/BRCA-interacting gene P-LP variants, and 3) No DNA repair germline variants.
  - Bar plots of GSVA score differences in HGG patients with versus without P-LP variants in DNA repair genes. 
  - Violin plots of GSVA scores in HGG patients with 1) MMR gene P-LP variants, 2) BRCA/BRCA-interacting gene P-LP variants, and 3) No DNA repair germline variants.
  
  
## Directory structure

```
.
├── 01-get-dna-repair-variant-samples.Rmd
├── 01-get-dna-repair-variant-samples.html
├── 02-mutational-signatures.Rmd
├── 02-mutational-signatures.html
├── 03-run-gsva-comparison.Rmd
├── 03-run-gsva-comparison.html
├── README.md
├── input
│   ├── COSMICv3.3_signature_exposures_SBS39_included.tsv
│   ├── Knijnenburg_paper
│   ├── broad
│   └── pbta-hallmark-gsva-scores.tsv
├── plots
│   ├── Broad_GO_repair_genes
│   ├── DIPG or DMG
│   ├── HGG, H3 wildtype
│   ├── Knijnenburg_repair_genes
│   ├── Mesenchymal tumor
│   └── Non-neoplastic tumor
├── results
│   ├── Broad_GO_repair_genes
│   ├── DIPG or DMG
│   ├── HGG, H3 wildtype
│   ├── Knijnenburg_repair_genes
│   ├── Mesenchymal tumor
│   ├── Non-neoplastic tumor
│   └── pbta-KEGG-gsva-scores.tsv
├── run-kegg-gsva.R
├── run_dna_repair_variant_analysis.sh
└── util
    ├── gsva_functions.R
    └── mutsigs_functions.R
```
  