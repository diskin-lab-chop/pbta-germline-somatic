library(tidyverse) 
library(qs)

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "methylation")
input_dir <- file.path(analysis_dir, "input")

# set file paths

methyl_file <- file.path(data_dir,
                         "IlluminaHumanMethylationEPIC-methyl-beta-values-masked.rds")

methyl_manifest <- read_tsv(file.path(input_dir,
                                      "itt-141-pbta-methylation.tsv")) %>%
  dplyr::mutate(name = str_remove_all(file_name, 
                                      "CBTN-Methylation/|_Red.idat.gz|_Grn.idat.gz|_Grn.idat|_Red.idat")) %>%
  distinct(name, `Kids First Biospecimen ID`)

# Load methylation data
methyl <- readRDS(methyl_file) %>%
  column_to_rownames("Probe_ID")

# convert column names to BS IDs by matching file names to IDs in manifest
match_ids <- as.vector(unlist(methyl_manifest[match(colnames(methyl), methyl_manifest$name), "Kids First Biospecimen ID"]))

colnames(methyl) <- match_ids

#confirm rename
methyl[1:5,1:5]
names(methyl)

qsave(methyl, file.path(data_dir, "methyl-beta-values-masked.qs"), preset = "fast")

