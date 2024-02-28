# Assess cancer predisposition syndrome-associated gene P/LP variant frequency in PBTA germline cohort

This module identifies all P/LP variants in genes associated with cancer predisposition syndromes, and cross-references findings with clinical reports of syndromes to determine if findings are incidental. Conversely, the frequency of patients with a reported syndrome AND that have an identified P/LP variant in the associated gene is plotted. 

## Usage

`bash run_module` 

## Folder content 

1. `01-incidental-findings.Rmd`; pulls all variants in syndrome-associated genes, and determines if finding is incidental by integrating histology data

2. `02-plot-predispositions.Rmd`; creates predisposition variant plots

3. `input/` directory contains the following files:
  - `predisposition-syndromes.tsv`; list of syndromes and associated genes
  - `predispositions-path-review.tsv`; path report notes in patients with incidental findings

3. `results/` directory contains the following files: 
  - `incidental-findings-predisposition-variants.tsv`; all syndrome-associated snv/indel variants and initial assessment of whether incidental 
  - `incidental-findings-predisposition-variants-path-review.tsv`; update of above file, with incidental status updated after path report review
  - `incidental-findings-predisposition-structural-variants.tsv`; all syndrome-associated SVs and initial assessment of whether incidental 
  - `incidental-findings-predisposition-structural-variants-path-review.tsv`; update of above file, with incidental status updated after path report review
  - `predisposition-patients-no-plp-variants.tsv`; patients with reported syndrome but no P/LP variant identified
  
4. `plots/` directory contains the following figures: 
  - Bar plot of number of patients with identified P/LP variant by syndrome (`gene-plp-prevalence-by-syndrome.pdf`)
  - Bar plot of number of patietns with reported syndrome by syndrome-associated gene with identified P/LP variant (`syndrome-prevalence-by-gene.pdf`)
