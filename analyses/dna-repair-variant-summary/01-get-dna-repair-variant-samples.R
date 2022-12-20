
library(tidyverse)
library(data.table)

setwd('OpenPBTA-germline')

```{r load libraries}
library(data.table)
library(tidyverse)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "dna-repair-variant-summary")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")

histologies_file <- file.path(root_dir, 'analyses', 'collapse-tumor-histologies', 'results',"germline-primary-plus-tumor-histologies-CBTN-dx.tsv")
plotgroup_file <- file.path(input_dir, "hist_with_plotGroup.csv")

plp_file <- file.path(data_dir, "pbta_germline_plp_calls.tsv")

dna_repair_file <- file.path(input_dir, 'dna_repair_all.txt')
ber_file <- file.path(input_dir, 'base_excision_repair.txt')
hr_file <- file.path(input_dir, 'homologous_recombination.txt')
mmr_file <- file.path(input_dir, 'mismatch_repair.txt')
ner_file <- file.path(input_dir, 'nucleotide_excision_repair.txt')
nhej_file <- file.path(input_dir, 'nonhomologous_end_joining.txt')

cpg_file <- file.path(root_dir, 'analyses', 'oncokb-annotation', 'input', 'cpg.txt')
```

repair_genes <- read_tsv(dna_repair_file, col_names = F) %>%
  rename('genes' = X1)
ber_genes <- read_tsv(ber_file, col_names = F) %>%
  rename('genes' = X1)
hr_genes <- read_tsv(hr_file, col_names = F) %>%
  rename('genes' = X1)
mmr_genes <- read_tsv(mmr_file, col_names = F) %>%
  rename('genes' = X1)
ner_genes <- read_tsv(ner_file, col_names = F) %>%
  rename('genes' = X1)
nhej_genes <- read_tsv(nhej_file, col_names = F) %>%
  rename('genes' = X1)

cpgs <- read_tsv(cpg_file, col_names = F) %>%
  rename('genes' = X1)

hist <- read_tsv(histologies_file)

germline <- read_tsv(plp_file) %>%
  filter(Kids_First_Biospecimen_ID %in% hist$Kids_First_Biospecimen_ID_normal)


repair_germline <- germline %>%
  filter(Hugo_Symbol %in% repair_genes$genes & Hugo_Symbol %in% cpgs$genes) %>%
  add_column(class = 'DNA_repair') %>%
  readr::write_tsv(file.path(results_dir, 'dna_repair_germline_variant_plp.tsv'))

ber_germline <- germline %>%
  filter(Hugo_Symbol %in% ber_genes$genes & Hugo_Symbol %in% cpgs$genes) %>%
  add_column(class = rep('BER', nrow(ber_germline))) %>%
  readr::write_tsv(file.path(results_dir, 'ber_germline_variant_plp.tsv'))

hr_germline <- germline %>%
  filter(Hugo_Symbol %in% hr_genes$genes & Hugo_Symbol %in% cpgs$genes) %>%
  add_column(class = rep('HR', nrow(hr_germline))) %>%
  readr::write_tsv(file.path(results_dir, 'hr_germline_variant_plp.tsv'))

mmr_germline <- germline %>%
  filter(Hugo_Symbol %in% mmr_genes$genes & Hugo_Symbol %in% cpgs$genes) %>%
  add_column(class = rep('MMR', nrow(mmr_germline))) %>%
  readr::write_tsv(file.path(results_dir, 'mmr_germline_variant_plp.tsv'))

ner_germline <- germline %>%
  filter(Hugo_Symbol %in% ner_genes$genes & Hugo_Symbol %in% cpgs$genes) %>%
  add_column(class = rep('NER', nrow(ner_germline))) %>%
  readr::write_tsv(file.path(results_dir, 'ner_germline_variant_plp.tsv'))

nhej_germline <- germline %>%
  filter(Hugo_Symbol %in% nhej_genes$genes & Hugo_Symbol %in% cpgs$genes) %>%
  add_column(class = rep('NHEJ', nrow(nhej_germline))) %>%
  readr::write_tsv(file.path(results_dir, 'nhej_germline_variant_plp.tsv'))


dna_repair_summary <- germline %>%
  select(Kids_First_Biospecimen_ID) %>%
  distinct() %>%
  mutate(
    dna_repair_germline = ifelse(Kids_First_Biospecimen_ID %in% repair_germline$Kids_First_Biospecimen_ID,
                                      'Yes', 'No'),
    ber_germline = ifelse(Kids_First_Biospecimen_ID %in% ber_germline$Kids_First_Biospecimen_ID,
                                 'Yes', 'No'),
    hr_germline = ifelse(Kids_First_Biospecimen_ID %in% hr_germline$Kids_First_Biospecimen_ID,
                                 'Yes', 'No'),
    mmr_germline = ifelse(Kids_First_Biospecimen_ID %in% mmr_germline$Kids_First_Biospecimen_ID,
                                 'Yes', 'No'),
    ner_germline = ifelse(Kids_First_Biospecimen_ID %in% ner_germline$Kids_First_Biospecimen_ID,
                                 'Yes', 'No'),
    nhej_germline = ifelse(Kids_First_Biospecimen_ID %in% nhej_germline$Kids_First_Biospecimen_ID,
                                 'Yes', 'No')
    ) %>%
  readr::write_tsv(file.path(results_dir, 'dna_repair_germline_plp_summary.tsv'))




plotgroup <- read_csv(plotgroup_file)
head(plotgroup)

res <- repair_germline %>%
  bind_rows(ber_germline) %>%
  bind_rows(hr_germline) %>%
  bind_rows(mmr_germline) %>%
  bind_rows(ner_germline) %>%
  bind_rows(nhej_germline) %>%
  rename(Kids_First_Biospecimen_ID_normal = Kids_First_Biospecimen_ID) %>%
  left_join(plotgroup[,c("Kids_First_Biospecimen_ID_normal", "plot_group", "broad_histology_hex")], by = 'Kids_First_Biospecimen_ID_normal') %>%
  count(class, plot_group, broad_histology_hex)

res <- repair_germline %>%
  count(class, plot_group)

res$class <- factor(res$class, levels = c('DNA_repair', 'BER', 'HR', 'MMR', 'NER', 'NHEJ'))
res$plot_group <- factor(res$plot_group, levels = unique(res$plot_group))
res$broad_histology_hex <- factor(res$broad_histology_hex, levels = unique(res$broad_histology_hex))

cols <- unique(res[!duplicated(res$plot_group),]$broad_histology_hex)
cols

ggplot(res, aes(x = class, y = n, fill = plot_group)) +
  geom_bar(position="stack", stat="identity") +
  labs(y = 'count', x = 'germline variant class') + 
  theme_minimal()
  scale_fill_manual(values = unique(res$broad_histology_hex))

