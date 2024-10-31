## Plot PBTA ancestry prediction from Somalier analysis

## Ryan Corbett
## 2023


#### Load packagers

library(tidyverse)
library(ggplot2)
library(ggpubr)

## Set directory paths

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "demo-clin-stats")
plots_dir <- file.path(analysis_dir, "plots")

# If the plots directory does not exist, create it
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Call plotting theme script
source(file.path(root_dir, "figures", "theme.R"))


## Set file paths
cbtn_histologies_file <- file.path(root_dir, "analyses", "collapse-tumor-histologies", "results", "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

## Read histologies file
hist <- read_tsv(cbtn_histologies_file)

# define colorblind-friendly palette
okabe_palette <- colorblindr::palette_OkabeIto[c(1:3,5:6)]

# assign predicted ancestries to each palette color: African (AFR), Admixed American (AMR), East Asian (EAS), European (EUR), and South Asian (SAS)
names(okabe_palette) <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Plot PC1 and PC2 from Somalier ancestry prediction
pdf(file.path(plots_dir, "predicted-ancestry-pca.pdf"),
    height = 4, width = 10)

pc12 <- hist %>%
  dplyr::filter(!is.na(predicted_ancestry)) %>%
  ggplot(aes(x = PC1, y = PC2, fill = predicted_ancestry)) +
  geom_point(size=2, shape=23,
             show.legend = FALSE) +
  scale_fill_manual(values = okabe_palette,
                    labels=c(glue::glue("AFR (n={sum(grepl('AFR', hist$predicted_ancestry))})"),
                             glue::glue("AMR (n={sum(grepl('AMR', hist$predicted_ancestry))})"),
                             glue::glue("EAS (n={sum(grepl('EAS', hist$predicted_ancestry))})"),
                             glue::glue("EUR (n={sum(grepl('EUR', hist$predicted_ancestry))})"),
                             glue::glue("SAS (n={sum(grepl('SAS', hist$predicted_ancestry))})"))) +
  labs(fill = "predictd ancestry") +
  theme_Publication()

# Plot PC3 and PC4
pc34 <- hist %>%
  dplyr::filter(!is.na(predicted_ancestry)) %>%
  ggplot(aes(x = PC3, y = PC4, fill = predicted_ancestry)) +
  geom_point(size=2, shape=23) +
  scale_fill_manual(values = okabe_palette,
                    labels=c(glue::glue("AFR (n={sum(grepl('AFR', hist$predicted_ancestry))})"),
                             glue::glue("AMR (n={sum(grepl('AMR', hist$predicted_ancestry))})"),
                             glue::glue("EAS (n={sum(grepl('EAS', hist$predicted_ancestry))})"),
                             glue::glue("EUR (n={sum(grepl('EUR', hist$predicted_ancestry))})"),
                             glue::glue("SAS (n={sum(grepl('SAS', hist$predicted_ancestry))})"))) +
  labs(fill = "predicted ancestry") +
  theme_Publication()

ggarrange(pc12, pc34,
          widths = c(1.3,2.0))

dev.off()
