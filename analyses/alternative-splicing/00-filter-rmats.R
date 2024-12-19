# Filter splicing rmats data for downstream alternative splicing analyses
#
# Ryan Corbett
#
# 2024

library(tidyverse) 
library(data.table)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "alternative-splicing")
results_dir <- file.path(analysis_dir, "results")

# define file paths

hist_file <- file.path(data_dir,
                  "histologies.tsv")

splice_events_file <- file.path(data_dir,
                                "splice-events-rmats.tsv.gz")

cpg_file <- file.path(root_dir, "analyses",
                  "oncokb-annotation", "input",
                  "cpg.txt")

plp_file <- file.path(data_dir, 
                      "pbta-merged-plp-variants-autogvp-abridged.tsv")

plp_lowVAF_file <- file.path(data_dir,
                             "pbta-merged-plp-variants-autogvp-abridged-lowVAF-cpgs.tsv")

plp_sv_file <- file.path(data_dir,
                         "pbta_germline_svs.tsv")

# wrangle data

cpgs <- read_lines(cpg_file)

plp <- read_tsv(plp_file) %>%
  bind_rows(read_tsv(plp_lowVAF_file)) %>%
  dplyr::filter(gene_symbol_vep %in% cpgs)

plp_sv <- read_tsv(plp_sv_file)

# Load rmats and filter for CPGs in plp files

splice_events_cpgs <- fread(splice_events_file) %>%
  dplyr::filter(geneSymbol %in% c(plp$gene_symbol_vep,
                                  plp_sv$Hugo_Symbol_cpg),
                splicing_case != "MXE")

# write output

write_tsv(splice_events_cpgs,
          file.path(results_dir, "splice-events-rmats-cpgs.tsv.gz"))
