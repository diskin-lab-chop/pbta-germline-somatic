
<!-- rnb-text-begin -->

---
title: "Calculate somatic hit enrichment among cancer predisposition genes"
output: 
  html_notebook:
  toc: TRUE
toc_float: TRUE
author: Ryan Corbett
---

This script tests for enrichment of somatic alterations in cancer predisposition genes among patients with germline P-LP variants in the same gene, relative to the rest of the cohort

# Load libraries

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbInN1cHByZXNzUGFja2FnZVN0YXJ0dXBNZXNzYWdlcyh7XG4gIGxpYnJhcnkoZ2dwbG90MilcbiAgbGlicmFyeSh0aWR5dmVyc2UpXG4gIGxpYnJhcnkoQ29tcGxleEhlYXRtYXApXG4gIGxpYnJhcnkoY2lyY2xpemUpXG4gIGxpYnJhcnkoZ2dwdWJyKVxufSkiLCIiLCIjIFNldCB1cCBkaXJlY3RvcmllcyIsInJvb3RfZGlyIDwtIHJwcm9qcm9vdDo6ZmluZF9yb290KHJwcm9qcm9vdDo6aGFzX2RpcihcIi5naXRcIikpIiwiIiwiZGF0YV9kaXIgPC0gZmlsZS5wYXRoKHJvb3RfZGlyLCBcImRhdGFcIikiLCJhbmFseXNpc19kaXIgPC0gZmlsZS5wYXRoKHJvb3RfZGlyLCBcImFuYWx5c2VzXCIsIFwidHdvLWhpdHNcIikiLCJwbG90X2RpciA8LSBmaWxlLnBhdGgoYW5hbHlzaXNfZGlyLCBcInBsb3RzXCIpIiwicmVzdWx0c19kaXIgPC0gZmlsZS5wYXRoKGFuYWx5c2lzX2RpciwgXCJyZXN1bHRzXCIpIiwiIiwiIyBJZiB0aGUgZGlyZWN0b3J5IGRvZXMgbm90IGV4aXN0LCBjcmVhdGUgaXQiLCJpZiAoIWRpci5leGlzdHMocmVzdWx0c19kaXIpKSB7XG4gIGRpci5jcmVhdGUocmVzdWx0c19kaXIsIHJlY3Vyc2l2ZSA9IFRSVUUpXG59IiwiIiwiaWYgKCFkaXIuZXhpc3RzKHBsb3RfZGlyKSkge1xuICBkaXIuY3JlYXRlKHBsb3RfZGlyLCByZWN1cnNpdmUgPSBUUlVFKVxufSIsIiIsIiMgc291cmNlIHB1YmxpY2F0aW9uIHRoZW1lIiwic291cmNlKGZpbGUucGF0aChyb290X2RpciwgXCJmaWd1cmVzXCIsIFwidGhlbWUuUlwiKSkiXX0= -->

```r
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
})

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "two-hits")
plot_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")

# If the directory does not exist, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# source publication theme
source(file.path(root_dir, "figures", "theme.R"))
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Set input file paths

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbInBscF9maWxlIDwtIGZpbGUucGF0aChkYXRhX2RpciwgXCJwYnRhX2dlcm1saW5lX3BscF9jYWxscy50c3ZcIikiLCIiLCJoaXN0X2ZpbGUgPC0gZmlsZS5wYXRoKGRhdGFfZGlyLCBcImhpc3RvbG9naWVzLnRzdlwiKSIsIiIsIm1hZl9maWxlIDwtIGZpbGUucGF0aChyb290X2RpciwgXCJhbmFseXNlc1wiLCBcIm9uY29rYi1hbm5vdGF0aW9uXCIsIFwicmVzdWx0c1wiLCBcInNudi1jb25zZW5zdXMtcGx1cy1ob3RzcG90cy1nb2ktb25jb2tiLm1hZi50c3ZcIikiLCJjbnZfZmlsZSA8LSBmaWxlLnBhdGgoZGF0YV9kaXIsIFwiY29uc2Vuc3VzX3NlZ19hbm5vdGF0ZWRfY25fYXV0b3NvbWVzLnRzdi5nelwiKSIsImxvaF9maWxlIDwtIGZpbGUucGF0aChkYXRhX2RpciwgXCJwYnRhLWFsbC1nZW5lLWxvaC50c3YuZ3pcIikiLCJleHByX2ZpbGUgPC0gZmlsZS5wYXRoKGRhdGFfZGlyLCBcImdlbmUtZXhwcmVzc2lvbi1yc2VtLXRwbS1jb2xsYXBzZWQucmRzXCIpIiwiIiwic3BlY2ltZW5zX2ZpbGUgPC0gZmlsZS5wYXRoKGRhdGFfZGlyLCBcImluZGVwZW5kZW50LXNwZWNpbWVucy5ybmFzZXFwYW5lbC5wcmltYXJ5LXBsdXMudHN2XCIpIl19 -->

```r
plp_file <- file.path(data_dir, "pbta_germline_plp_calls.tsv")

hist_file <- file.path(data_dir, "histologies.tsv")

maf_file <- file.path(root_dir, "analyses", "oncokb-annotation", "results", "snv-consensus-plus-hotspots-goi-oncokb.maf.tsv")
cnv_file <- file.path(data_dir, "consensus_seg_annotated_cn_autosomes.tsv.gz")
loh_file <- file.path(data_dir, "pbta-all-gene-loh.tsv.gz")
expr_file <- file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds")

specimens_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary-plus.tsv")
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


# Load files

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbIiMgcmVhZCBpbiBoaXN0b2xvZ3kgZmlsZSIsImhpc3RfY29ob3J0IDwtIHJlYWRfdHN2KGZpbGUucGF0aChyb290X2RpciwgXCJhbmFseXNlc1wiLCBcImNvbGxhcHNlLXR1bW9yLWhpc3RvbG9naWVzXCIsIFwicmVzdWx0c1wiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJnZXJtbGluZS1wcmltYXJ5LXBsdXMtdHVtb3ItaGlzdG9sb2dpZXMtcGxvdC1ncm91cHMtY2xpbi1tZXRhLnRzdlwiKSwgc2hvd19jb2xfdHlwZXMgPSBGQUxTRSkiLCIiLCJoaXN0X2FsbCA8LSByZWFkX3RzdihoaXN0X2ZpbGUsIGd1ZXNzX21heCA9IDEwMDAwMDApIl19 -->

```r
# read in histology file
hist_cohort <- read_tsv(file.path(root_dir, "analyses", "collapse-tumor-histologies", "results",
                           "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv"), show_col_types = FALSE)

hist_all <- read_tsv(hist_file, guess_max = 1000000)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiUm93czogNDM3MjkgQ29sdW1uczogNTVcbuKUgOKUgCBDb2x1bW4gc3BlY2lmaWNhdGlvbiDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIBcbkRlbGltaXRlcjogXCJcXHRcIlxuY2hyICgzOSk6IEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSUQsIHNhbXBsZV9pZCwgYWxpcXVvdF9pZCwgS2lkc19GaXJzdF9QYXJ0aS4uLlxuZGJsICgxNik6IGFnZV9hdF9kaWFnbm9zaXNfZGF5cywgT1NfZGF5cywgRUZTX2RheXMsIGFnZV9sYXN0X3VwZGF0ZV9kYXlzLCBuby4uLlxuXG7ihLkgVXNlIGBzcGVjKClgIHRvIHJldHJpZXZlIHRoZSBmdWxsIGNvbHVtbiBzcGVjaWZpY2F0aW9uIGZvciB0aGlzIGRhdGEuXG7ihLkgU3BlY2lmeSB0aGUgY29sdW1uIHR5cGVzIG9yIHNldCBgc2hvd19jb2xfdHlwZXMgPSBGQUxTRWAgdG8gcXVpZXQgdGhpcyBtZXNzYWdlLlxuIn0= -->

```
Rows: 43729 Columns: 55
── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (39): Kids_First_Biospecimen_ID, sample_id, aliquot_id, Kids_First_Parti...
dbl (16): age_at_diagnosis_days, OS_days, EFS_days, age_last_update_days, no...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjpbImhpc3RfYWxsX2dlbm9taWMgPC0gaGlzdF9hbGwgJT4lXG4gIGRwbHlyOjpmaWx0ZXIoc2FtcGxlX3R5cGUgPT0gXCJUdW1vclwiICYgZXhwZXJpbWVudGFsX3N0cmF0ZWd5ICVpbiUgYyhcIldHU1wiLCBcIldYU1wiKSkiLCIiLCJoaXN0X2FsbF9ybmEgPC0gaGlzdF9hbGwgJT4lXG4gIGRwbHlyOjpmaWx0ZXIoc2FtcGxlX3R5cGUgPT0gXCJUdW1vclwiICYgZXhwZXJpbWVudGFsX3N0cmF0ZWd5ID09IFwiUk5BLVNlcVwiKSIsIiIsInR1bW9yX25vcm1hbF9nZW5vbWVfZGYgPC0gaGlzdF9jb2hvcnQgJT4lXG4gIGRwbHlyOjpzZWxlY3QoLUtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfdHVtb3IpICU+JVxuICBsZWZ0X2pvaW4oaGlzdF9hbGxfZ2Vub21pY1ssYyhcIktpZHNfRmlyc3RfUGFydGljaXBhbnRfSURcIiwgXCJLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEXCIpXSkgJT4lXG4gIGRwbHlyOjpyZW5hbWUoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF90dW1vciA9IEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSUQpICU+JVxuICBkcGx5cjo6c2VsZWN0KEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsLCBLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX3R1bW9yKSJdfQ== -->

```r
hist_all_genomic <- hist_all %>%
  dplyr::filter(sample_type == "Tumor" & experimental_strategy %in% c("WGS", "WXS"))

hist_all_rna <- hist_all %>%
  dplyr::filter(sample_type == "Tumor" & experimental_strategy == "RNA-Seq")

tumor_normal_genome_df <- hist_cohort %>%
  dplyr::select(-Kids_First_Biospecimen_ID_tumor) %>%
  left_join(hist_all_genomic[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID")]) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Biospecimen_ID_normal, Kids_First_Biospecimen_ID_tumor)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiSm9pbmluZyB3aXRoIGBieSA9IGpvaW5fYnkoS2lkc19GaXJzdF9QYXJ0aWNpcGFudF9JRClgXG4ifQ== -->

```
Joining with `by = join_by(Kids_First_Participant_ID)`
```



<!-- rnb-message-end -->

<!-- rnb-warning-begin eyJkYXRhIjoiV2FybmluZyBpbiBsZWZ0X2pvaW4oLiwgaGlzdF9hbGxfZ2Vub21pY1ssIGMoXCJLaWRzX0ZpcnN0X1BhcnRpY2lwYW50X0lEXCIsIDogRWFjaCByb3cgaW4gYHhgIGlzIGV4cGVjdGVkIHRvIG1hdGNoIGF0IG1vc3QgMSByb3cgaW4gYHlgLlxu4oS5IFJvdyAxIG9mIGB4YCBtYXRjaGVzIG11bHRpcGxlIHJvd3MuXG7ihLkgSWYgbXVsdGlwbGUgbWF0Y2hlcyBhcmUgZXhwZWN0ZWQsIHNldCBgbXVsdGlwbGUgPSBcImFsbFwiYCB0byBzaWxlbmNlIHRoaXMgd2FybmluZy5cbiJ9 -->

```
Warning in left_join(., hist_all_genomic[, c("Kids_First_Participant_ID", : Each row in `x` is expected to match at most 1 row in `y`.
ℹ Row 1 of `x` matches multiple rows.
ℹ If multiple matches are expected, set `multiple = "all"` to silence this warning.
```



<!-- rnb-warning-end -->

<!-- rnb-source-begin eyJkYXRhIjoidHVtb3Jfbm9ybWFsX3JuYV9kZiA8LSBoaXN0X2NvaG9ydCAlPiVcbiAgbGVmdF9qb2luKGhpc3RfYWxsX3JuYVssYyhcIktpZHNfRmlyc3RfUGFydGljaXBhbnRfSURcIiwgXCJLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEXCIpXSkgJT4lXG4gIGRwbHlyOjpyZW5hbWUoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ybmEgPSBLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEKSAlPiVcbiAgZHBseXI6OnNlbGVjdChLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCwgS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ybmEpIn0= -->

```r
tumor_normal_rna_df <- hist_cohort %>%
  left_join(hist_all_rna[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID")]) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_rna = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Biospecimen_ID_normal, Kids_First_Biospecimen_ID_rna)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiSm9pbmluZyB3aXRoIGBieSA9IGpvaW5fYnkoS2lkc19GaXJzdF9QYXJ0aWNpcGFudF9JRClgXG4ifQ== -->

```
Joining with `by = join_by(Kids_First_Participant_ID)`
```



<!-- rnb-message-end -->

<!-- rnb-warning-begin eyJkYXRhIjoiV2FybmluZyBpbiBsZWZ0X2pvaW4oLiwgaGlzdF9hbGxfcm5hWywgYyhcIktpZHNfRmlyc3RfUGFydGljaXBhbnRfSURcIiwgXCJLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEXCIpXSk6IEVhY2ggcm93IGluIGB4YCBpcyBleHBlY3RlZCB0byBtYXRjaCBhdCBtb3N0IDEgcm93IGluIGB5YC5cbuKEuSBSb3cgMSBvZiBgeGAgbWF0Y2hlcyBtdWx0aXBsZSByb3dzLlxu4oS5IElmIG11bHRpcGxlIG1hdGNoZXMgYXJlIGV4cGVjdGVkLCBzZXQgYG11bHRpcGxlID0gXCJhbGxcImAgdG8gc2lsZW5jZSB0aGlzIHdhcm5pbmcuXG4ifQ== -->

```
Warning in left_join(., hist_all_rna[, c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID")]): Each row in `x` is expected to match at most 1 row in `y`.
ℹ Row 1 of `x` matches multiple rows.
ℹ If multiple matches are expected, set `multiple = "all"` to silence this warning.
```



<!-- rnb-warning-end -->

<!-- rnb-source-begin eyJkYXRhIjoibWFmIDwtIHJlYWRfdHN2KG1hZl9maWxlKSJ9 -->

```r
maf <- read_tsv(maf_file)
```



<!-- rnb-source-end -->

<!-- rnb-warning-begin eyJkYXRhIjoiV2FybmluZzogT25lIG9yIG1vcmUgcGFyc2luZyBpc3N1ZXMsIGNhbGwgYHByb2JsZW1zKClgIG9uIHlvdXIgZGF0YSBmcmFtZSBmb3IgZGV0YWlscywgZS5nLjpcbiAgZGF0IDwtIHZyb29tKC4uLilcbiAgcHJvYmxlbXMoZGF0KVxuIn0= -->

```
Warning: One or more parsing issues, call `problems()` on your data frame for details, e.g.:
  dat <- vroom(...)
  problems(dat)
```



<!-- rnb-warning-end -->

<!-- rnb-message-begin eyJkYXRhIjoiUm93czogNTk4NjggQ29sdW1uczogMTQzXG7ilIDilIAgQ29sdW1uIHNwZWNpZmljYXRpb24g4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSAXG5EZWxpbWl0ZXI6IFwiXFx0XCJcbmNociAoNzkpOiBIdWdvX1N5bWJvbCwgQ2VudGVyLCBOQ0JJX0J1aWxkLCBDaHJvbW9zb21lLCBTdHJhbmQsIFZhcmlhbnRfQ2xhc3MuLi5cbmRibCAoMzcpOiBFbnRyZXpfR2VuZV9JZCwgU3RhcnRfUG9zaXRpb24sIEVuZF9Qb3NpdGlvbiwgdF9kZXB0aCwgdF9yZWZfY291bnQuLi5cbm51bSAgKDEpOiBQVUJNRURcbmxnbCAoMjYpOiBkYlNOUF9WYWxfU3RhdHVzLCBUdW1vcl9WYWxpZGF0aW9uX0FsbGVsZTEsIFR1bW9yX1ZhbGlkYXRpb25fQWxsZWwuLi5cblxu4oS5IFVzZSBgc3BlYygpYCB0byByZXRyaWV2ZSB0aGUgZnVsbCBjb2x1bW4gc3BlY2lmaWNhdGlvbiBmb3IgdGhpcyBkYXRhLlxu4oS5IFNwZWNpZnkgdGhlIGNvbHVtbiB0eXBlcyBvciBzZXQgYHNob3dfY29sX3R5cGVzID0gRkFMU0VgIHRvIHF1aWV0IHRoaXMgbWVzc2FnZS5cbiJ9 -->

```
Rows: 59868 Columns: 143
── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (79): Hugo_Symbol, Center, NCBI_Build, Chromosome, Strand, Variant_Class...
dbl (37): Entrez_Gene_Id, Start_Position, End_Position, t_depth, t_ref_count...
num  (1): PUBMED
lgl (26): dbSNP_Val_Status, Tumor_Validation_Allele1, Tumor_Validation_Allel...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjoicGxwIDwtIHJlYWRfdHN2KHBscF9maWxlKSAlPiVcbiAgZHBseXI6OmZpbHRlcihLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEICVpbiUgaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwpIn0= -->

```r
plp <- read_tsv(plp_file) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% hist_cohort$Kids_First_Biospecimen_ID_normal)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiUm93czogNTY1OCBDb2x1bW5zOiA3XG7ilIDilIAgQ29sdW1uIHNwZWNpZmljYXRpb24g4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSAXG5EZWxpbWl0ZXI6IFwiXFx0XCJcbmNociAoNyk6IEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSUQsIGV4cGVyaW1lbnRhbF9zdHJhdGVneSwgTC1MUF9DYWxsX2ZpbmFsLCAuLi5cblxu4oS5IFVzZSBgc3BlYygpYCB0byByZXRyaWV2ZSB0aGUgZnVsbCBjb2x1bW4gc3BlY2lmaWNhdGlvbiBmb3IgdGhpcyBkYXRhLlxu4oS5IFNwZWNpZnkgdGhlIGNvbHVtbiB0eXBlcyBvciBzZXQgYHNob3dfY29sX3R5cGVzID0gRkFMU0VgIHRvIHF1aWV0IHRoaXMgbWVzc2FnZS5cbiJ9 -->

```
Rows: 5658 Columns: 7
── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (7): Kids_First_Biospecimen_ID, experimental_strategy, L-LP_Call_final, ...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjoiY252IDwtIHJlYWRfdHN2KGNudl9maWxlKSAlPiVcbiAgZHBseXI6OmZpbHRlcihiaW9zcGVjaW1lbl9pZCAlaW4lIHR1bW9yX25vcm1hbF9nZW5vbWVfZGYkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF90dW1vcikgJT4lXG4gIGRwbHlyOjpyZW5hbWUoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF90dW1vciA9IGJpb3NwZWNpbWVuX2lkKSAlPiVcbiAgbGVmdF9qb2luKHR1bW9yX25vcm1hbF9nZW5vbWVfZGYpICU+JVxuICBkaXN0aW5jdCgpIn0= -->

```r
cnv <- read_tsv(cnv_file) %>%
  dplyr::filter(biospecimen_id %in% tumor_normal_genome_df$Kids_First_Biospecimen_ID_tumor) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = biospecimen_id) %>%
  left_join(tumor_normal_genome_df) %>%
  distinct()
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiTXVsdGlwbGUgZmlsZXMgaW4gemlwOiByZWFkaW5nICdjb25zZW5zdXNfc2VnX2Fubm90YXRlZF9jbl9hdXRvc29tZXMudHN2J1xuUm93czogMzkxNzAwMCBDb2x1bW5zOiA3XG7ilIDilIAgQ29sdW1uIHNwZWNpZmljYXRpb24g4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSAXG5EZWxpbWl0ZXI6IFwiXFx0XCJcbmNociAoNSk6IGJpb3NwZWNpbWVuX2lkLCBzdGF0dXMsIGVuc2VtYmwsIGdlbmVfc3ltYm9sLCBjeXRvYmFuZFxuZGJsICgyKTogY29weV9udW1iZXIsIHBsb2lkeVxuXG7ihLkgVXNlIGBzcGVjKClgIHRvIHJldHJpZXZlIHRoZSBmdWxsIGNvbHVtbiBzcGVjaWZpY2F0aW9uIGZvciB0aGlzIGRhdGEuXG7ihLkgU3BlY2lmeSB0aGUgY29sdW1uIHR5cGVzIG9yIHNldCBgc2hvd19jb2xfdHlwZXMgPSBGQUxTRWAgdG8gcXVpZXQgdGhpcyBtZXNzYWdlLlxuSm9pbmluZyB3aXRoIGBieSA9IGpvaW5fYnkoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF90dW1vcilgXG4ifQ== -->

```
Multiple files in zip: reading 'consensus_seg_annotated_cn_autosomes.tsv'
Rows: 3917000 Columns: 7
── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (5): biospecimen_id, status, ensembl, gene_symbol, cytoband
dbl (2): copy_number, ploidy

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
Joining with `by = join_by(Kids_First_Biospecimen_ID_tumor)`
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjoibG9oIDwtIHJlYWRfdHN2KGxvaF9maWxlKSAlPiVcbiAgZHBseXI6OmZpbHRlcihLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCAlaW4lIHR1bW9yX25vcm1hbF9nZW5vbWVfZGYkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwpIn0= -->

```r
loh <- read_tsv(loh_file) %>%
  dplyr::filter(Kids_First_Biospecimen_ID_normal %in% tumor_normal_genome_df$Kids_First_Biospecimen_ID_normal)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiTXVsdGlwbGUgZmlsZXMgaW4gemlwOiByZWFkaW5nICdwYnRhLWFsbC1nZW5lLWxvaC50c3YnXG5Sb3dzOiAxNDk1MjQyNSBDb2x1bW5zOiA0XG7ilIDilIAgQ29sdW1uIHNwZWNpZmljYXRpb24g4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSAXG5EZWxpbWl0ZXI6IFwiXFx0XCJcbmNociAoMyk6IEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsLCBLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX3R1bW9yLCAuLi5cbmRibCAoMSk6IExPSF9zY29yZVxuXG7ihLkgVXNlIGBzcGVjKClgIHRvIHJldHJpZXZlIHRoZSBmdWxsIGNvbHVtbiBzcGVjaWZpY2F0aW9uIGZvciB0aGlzIGRhdGEuXG7ihLkgU3BlY2lmeSB0aGUgY29sdW1uIHR5cGVzIG9yIHNldCBgc2hvd19jb2xfdHlwZXMgPSBGQUxTRWAgdG8gcXVpZXQgdGhpcyBtZXNzYWdlLlxuIn0= -->

```
Multiple files in zip: reading 'pbta-all-gene-loh.tsv'
Rows: 14952425 Columns: 4
── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (3): Kids_First_Biospecimen_ID_normal, Kids_First_Biospecimen_ID_tumor, ...
dbl (1): LOH_score

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjpbImV4cHIgPC0gcmVhZFJEUyhleHByX2ZpbGUpIiwiIyBzdWJzZXQgYGV4cHJgIHRvIG9ubHkgaW5jbHVkZSBCUyBJRHMgaW4gaGlzdCIsInN1YnNldF9leHByIDwtIGV4cHJbLGludGVyc2VjdCh0dW1vcl9ub3JtYWxfcm5hX2RmJEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfcm5hLCBjb2xuYW1lcyhleHByKSldIl19 -->

```r
expr <- readRDS(expr_file)
# subset `expr` to only include BS IDs in hist
subset_expr <- expr[,intersect(tumor_normal_rna_df$Kids_First_Biospecimen_ID_rna, colnames(expr))]
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


calculate z-scores from TPM data

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbIiMgbG9nMih4ICsgMSkgdHJhbnNmb3JtIHRoZSBleHByZXNzaW9uIG1hdHJpeCIsImxvZ19leHByZXNzaW9uIDwtIGxvZzIoc3Vic2V0X2V4cHIgKyAxKSIsIiIsIiMgU2NhbGUgdGhlIGdlbmUgdmFsdWVzIC0tIHNjYWxlKCkgd29ya3Mgb24gdGhlIGNvbHVtbnMsIGhlbmNlIHRoZSB0cmFuc3Bvc2UiLCJ6X3Njb3JlZF9leHByZXNzaW9uIDwtIHNjYWxlKHQobG9nX2V4cHJlc3Npb24pLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjZW50ZXIgPSBUUlVFLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICBzY2FsZSA9IFRSVUUpIiwiel9zY29yZWRfZXhwcmVzc2lvbiA8LSBhcy5kYXRhLmZyYW1lKHpfc2NvcmVkX2V4cHJlc3Npb24pIiwiIiwiel9zY29yZWRfZXhwcmVzc2lvbiA8LSB6X3Njb3JlZF9leHByZXNzaW9uICU+JVxuICB0KCkgJT4lXG4gIGFzLmRhdGEuZnJhbWUoKSAlPiVcbiAgcm93bmFtZXNfdG9fY29sdW1uKFwiSHVnb19TeW1ib2xcIikgJT4lXG4gIGdhdGhlcihrZXkgPSBcIktpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURcIiwgdmFsdWUgPSBcImV4cHJfenNjb3JlXCIsIC1IdWdvX1N5bWJvbCkgJT4lXG4gIGRwbHlyOjpyZW5hbWUoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ybmEgPSBLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEKSAlPiVcbiAgbGVmdF9qb2luKHR1bW9yX25vcm1hbF9ybmFfZGYpIl19 -->

```r
# log2(x + 1) transform the expression matrix
log_expression <- log2(subset_expr + 1)

# Scale the gene values -- scale() works on the columns, hence the transpose
z_scored_expression <- scale(t(log_expression),
                             center = TRUE,
                             scale = TRUE)
z_scored_expression <- as.data.frame(z_scored_expression)

z_scored_expression <- z_scored_expression %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Hugo_Symbol") %>%
  gather(key = "Kids_First_Biospecimen_ID", value = "expr_zscore", -Hugo_Symbol) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_rna = Kids_First_Biospecimen_ID) %>%
  left_join(tumor_normal_rna_df)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiSm9pbmluZyB3aXRoIGBieSA9IGpvaW5fYnkoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ybmEpYFxuIn0= -->

```
Joining with `by = join_by(Kids_First_Biospecimen_ID_rna)`
```



<!-- rnb-message-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Identify cancer predisposition genes with P-LP variants in >= 5 patients

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbImNwZyA8LSByZWFkX2xpbmVzKGZpbGUucGF0aChyb290X2RpciwgXCJhbmFseXNlc1wiLCBcIm9uY29rYi1hbm5vdGF0aW9uXCIsIFwiaW5wdXRcIiwgXCJjcGcudHh0XCIpKSIsIiIsImdlbmVzIDwtIHBscCAlPiVcbiAgZHBseXI6OmNvdW50KEh1Z29fU3ltYm9sKSAlPiVcbiAgZHBseXI6OmZpbHRlcihuID49IDQgJiBIdWdvX1N5bWJvbCAlaW4lIGNwZykgJT4lXG4gIHB1bGwoSHVnb19TeW1ib2wpIl19 -->

```r
cpg <- read_lines(file.path(root_dir, "analyses", "oncokb-annotation", "input", "cpg.txt"))

genes <- plp %>%
  dplyr::count(Hugo_Symbol) %>%
  dplyr::filter(n >= 4 & Hugo_Symbol %in% cpg) %>%
  pull(Hugo_Symbol)
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Create empty data frames to store enrichment stats

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbImVtcHR5X2RmIDwtIGRhdGEuZnJhbWUocm93Lm5hbWVzID0gZ2VuZXMsXG4gICAgICAgICAgICAgICAgICAgICBhbHRlcmF0aW9uID0gcmVwKFwiIFwiLCBsZW5ndGgoZ2VuZXMpKSxcbiAgICAgICAgICAgICAgICAgICAgIG5fcGxwX2FsdCA9IHJlcCgwLCBsZW5ndGgoZ2VuZXMpKSwgXG4gICAgICAgICAgICAgICAgICAgICBuX3BscF9ub19hbHQgPSByZXAoMCwgbGVuZ3RoKGdlbmVzKSksXG4gICAgICAgICAgICAgICAgICAgICBuX25vX3BscF9hbHQgPSByZXAoMCwgbGVuZ3RoKGdlbmVzKSksXG4gICAgICAgICAgICAgICAgICAgICBuX25vX3BscF9ub19hbHQgPSByZXAoMCwgbGVuZ3RoKGdlbmVzKSksXG4gICAgICAgICAgICAgICAgICAgICBvZGRzUmF0aW8gPSByZXAoMCwgbGVuZ3RoKGdlbmVzKSksXG4gICAgICAgICAgICAgICAgICAgICBwID0gcmVwKDAsIGxlbmd0aChnZW5lcykpLFxuICAgICAgICAgICAgICAgICAgICAgQ0lfbG93ZXIgPSByZXAoMCwgbGVuZ3RoKGdlbmVzKSksXG4gICAgICAgICAgICAgICAgICAgICBDSV91cHBlciA9IHJlcCgwLCBsZW5ndGgoZ2VuZXMpKSxcbiAgICAgICAgICAgICAgICAgICAgIGFkanVzdGVkX3AgPSByZXAoMCwgbGVuZ3RoKGdlbmVzKSkpIiwiIiwic252X2RmIDwtIGVtcHR5X2RmICU+JVxuICBkcGx5cjo6bXV0YXRlKGFsdGVyYXRpb24gPSBcIk9uY29nZW5pYyBTTlZcIikiLCIiLCJjbnZfbG9zc19kZiA8LSBlbXB0eV9kZiAlPiVcbiAgbXV0YXRlKGFsdGVyYXRpb24gPSBcIkNOIGxvc3NcIikiLCIiLCJjbnZfZ2Fpbl9kZiA8LSBlbXB0eV9kZiAlPiVcbiAgZHBseXI6Om11dGF0ZShhbHRlcmF0aW9uID0gXCJDTiBnYWluXCIpIiwiIiwibG9oX2RmIDwtIGVtcHR5X2RmICU+JVxuICBkcGx5cjo6bXV0YXRlKGFsdGVyYXRpb24gPSBcIkxPSCA+IDAuMjVcIikiLCIiLCJleHByX2dhaW5fZGYgPC0gZW1wdHlfZGYgJT4lXG4gIGRwbHlyOjptdXRhdGUoYWx0ZXJhdGlvbiA9IFwiRXhwciBUUE0gPiAyXCIpIiwiIiwiZXhwcl9sb3NzX2RmIDwtIGVtcHR5X2RmICU+JVxuICBkcGx5cjo6bXV0YXRlKGFsdGVyYXRpb24gPSBcIkV4cHIgVFBNIDwgLTJcIikiXX0= -->

```r
empty_df <- data.frame(row.names = genes,
                     alteration = rep(" ", length(genes)),
                     n_plp_alt = rep(0, length(genes)), 
                     n_plp_no_alt = rep(0, length(genes)),
                     n_no_plp_alt = rep(0, length(genes)),
                     n_no_plp_no_alt = rep(0, length(genes)),
                     oddsRatio = rep(0, length(genes)),
                     p = rep(0, length(genes)),
                     CI_lower = rep(0, length(genes)),
                     CI_upper = rep(0, length(genes)),
                     adjusted_p = rep(0, length(genes)))

snv_df <- empty_df %>%
  dplyr::mutate(alteration = "Oncogenic SNV")

cnv_loss_df <- empty_df %>%
  mutate(alteration = "CN loss")

cnv_gain_df <- empty_df %>%
  dplyr::mutate(alteration = "CN gain")

loh_df <- empty_df %>%
  dplyr::mutate(alteration = "LOH > 0.25")

expr_gain_df <- empty_df %>%
  dplyr::mutate(alteration = "Expr TPM > 2")

expr_loss_df <- empty_df %>%
  dplyr::mutate(alteration = "Expr TPM < -2")
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Loop through genes and calculate enrichment of somatic alterations among patients with germline P-LP variant in same gene

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiZm9yIChnZW5lIGluIGdlbmVzKXtcbiAgXG4gICMgaWRlbnRpZnkgc2FtcGxlcyB3aXRoIFAtTFAgdmFyaWFudCBpbiBnZW5lXG4gIHBscF9nZW5lX2JzIDwtIHBscCAlPiVcbiAgICBkcGx5cjo6ZmlsdGVyKEh1Z29fU3ltYm9sID09IGdlbmUpICU+JVxuICAgIHB1bGwoS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRCkgXG5cbiAgIyBGaWx0ZXIgTUFGIGZvciBvbmx5IHZhcmlhbnRzIGluIGdlbmVcbiAgZ2VuZV9tYWYgPC0gbWFmICU+JVxuICAgICAgZHBseXI6OmZpbHRlcihIdWdvX1N5bWJvbCA9PSBnZW5lICYgT05DT0dFTklDICVpbiUgYyhcIk9uY29nZW5pY1wiLCBcIkxpa2VseSBPbmNvZ2VuaWNcIikpIFxuICBcbiAgIyBjcmVhdGUgY29udGluZ2VuY3kgdGFibGUgb2YgcGxwIFQvRiBhbmQgc29tYXRpYyBhbHRlcmF0aW9uIFQvRlxuICBzbnZfbWF0IDwtIHRhYmxlKGlzX3BscCA9IGhpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsICVpbiUgcGxwX2dlbmVfYnMsIFxuICAgICAgICAgICAgICAgICAgIGlzX2FsdCA9IGZhY3RvcihoaXN0X2NvaG9ydCRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCAlaW4lIGdlbmVfbWFmJE1hdGNoZWRfTm9ybV9TYW1wbGVfQmFyY29kZSwgbGV2ZWxzID0gYyhcIkZBTFNFXCIsIFwiVFJVRVwiKSkpXG4gIFxuICAjIGFkZCBzdGF0cyB0byBkZiBcbiAgc252X2RmW2dlbmUsIGMoXCJuX25vX3BscF9ub19hbHRcIiwgXCJuX3BscF9ub19hbHRcIiwgXCJuX25vX3BscF9hbHRcIiwgXCJuX3BscF9hbHRcIildID0gdW5saXN0KHNudl9tYXQpXG4gIHNudl9kZltnZW5lLCBjKFwib2Rkc1JhdGlvXCIsIFwicFwiLCBcIkNJX2xvd2VyXCIsIFwiQ0lfdXBwZXJcIildIDwtIHVubGlzdChmaXNoZXIudGVzdChzbnZfbWF0KVtjKFwiZXN0aW1hdGVcIiwgXCJwLnZhbHVlXCIsIFwiY29uZi5pbnRcIildKVxuICBcbiAgIyBmaWx0ZXIgY252IGRmIGZvciBzYW1wbGVzIHdpdGggbG9zcyBpbiBnZW5lXG4gIGdlbmVfY252X2xvc3MgPC0gY252ICU+JVxuICAgIGRwbHlyOjpmaWx0ZXIoZ2VuZV9zeW1ib2wgPT0gZ2VuZSAmIHN0YXR1cyA9PSBcImxvc3NcIilcbiAgXG4gICMgY3JlYXRlIFAvTFAgYW5kIHNvbWF0aWMgYWx0ZXJhdGlvbiBjb250aW5nZW5jeSB0YWJsZSBhbmQgc3RvcmUgc3RhdHMgaW4gZGZcbiAgY252X2xvc3NfbWF0IDwtIHRhYmxlKGlzX3BscCA9IGhpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsW2hpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsICVpbiUgY252JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsXSAlaW4lIHBscF9nZW5lX2JzLFxuICAgICAgICAgICAgICAgICAgICAgICAgaXNfYWx0ID0gZmFjdG9yKGhpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsW2hpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsICVpbiUgY252JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsXSAlaW4lIGdlbmVfY252X2xvc3MkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwsIGxldmVscyA9IGMoXCJGQUxTRVwiLCBcIlRSVUVcIikpKVxuICBcbiAgY252X2xvc3NfZGZbZ2VuZSwgYyhcIm5fbm9fcGxwX25vX2FsdFwiLCBcIm5fcGxwX25vX2FsdFwiLCBcIm5fbm9fcGxwX2FsdFwiLCBcIm5fcGxwX2FsdFwiKV0gPSB1bmxpc3QoY252X2xvc3NfbWF0KVxuICBjbnZfbG9zc19kZltnZW5lLCBjKFwib2Rkc1JhdGlvXCIsIFwicFwiLCBcIkNJX2xvd2VyXCIsIFwiQ0lfdXBwZXJcIildIDwtIHVubGlzdChmaXNoZXIudGVzdChjbnZfbG9zc19tYXQpW2MoXCJlc3RpbWF0ZVwiLCBcInAudmFsdWVcIiwgXCJjb25mLmludFwiKV0pXG4gIFxuICAjIFJlcGVhdCBmb3IgY29weSBudW1iZXIgZ2FpbnNcbiAgZ2VuZV9jbnZfZ2FpbiA8LSBjbnYgJT4lXG4gICAgZHBseXI6OmZpbHRlcihnZW5lX3N5bWJvbCA9PSBnZW5lICYgc3RhdHVzID09IFwiZ2FpblwiKVxuICBcbiAgY252X2dhaW5fbWF0IDwtIHRhYmxlKGlzX3BscCA9IGhpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsW2hpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsICVpbiUgY252JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsXSAlaW4lIHBscF9nZW5lX2JzLFxuICAgICAgICAgICAgICAgICAgICAgICAgaXNfYWx0ID0gZmFjdG9yKGhpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsW2hpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsICVpbiUgY252JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsXSAlaW4lIGdlbmVfY252X2dhaW4kS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwsIGxldmVscyA9IGMoXCJGQUxTRVwiLCBcIlRSVUVcIikpKVxuICBcbiAgY252X2dhaW5fZGZbZ2VuZSwgYyhcIm5fbm9fcGxwX25vX2FsdFwiLCBcIm5fcGxwX25vX2FsdFwiLCBcIm5fbm9fcGxwX2FsdFwiLCBcIm5fcGxwX2FsdFwiKV0gPSB1bmxpc3QoY252X2dhaW5fbWF0KVxuICBjbnZfZ2Fpbl9kZltnZW5lLCBjKFwib2Rkc1JhdGlvXCIsIFwicFwiLCBcIkNJX2xvd2VyXCIsIFwiQ0lfdXBwZXJcIildIDwtIHVubGlzdChmaXNoZXIudGVzdChjbnZfZ2Fpbl9tYXQpW2MoXCJlc3RpbWF0ZVwiLCBcInAudmFsdWVcIiwgXCJjb25mLmludFwiKV0pXG4gIFxuICAjIEZpbHRlciBMT0ggZGYgZm9yIHNhbXBsZXMgZXhoaWJpdGluZyBMT0ggaW4gZ2VuZVxuICBnZW5lX2xvaCA8LSBsb2ggJT4lXG4gICAgZHBseXI6OmZpbHRlcihIdWdvX1N5bWJvbCA9PSBnZW5lICYgTE9IX3Njb3JlID4gMC4yNSlcbiAgXG4gICMgQ3JlYXRlIGNvbnRpbmdlbmN5IHRhYmxlIGFuZCBzdG9yZSBzdGF0cyBmb3IgTE9IIGVucmljaG1lbnQgXG4gIGxvaF9tYXQgPC0gdGFibGUoaXNfcGxwID0gaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxbaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwgJWluJSBsb2gkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxdICVpbiUgcGxwX2dlbmVfYnMsXG4gICAgICAgICAgICAgICAgICAgICAgICBpc19hbHQgPSBmYWN0b3IoaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxbaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwgJWluJSBsb2gkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxdICVpbiUgZ2VuZV9sb2gkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwsIGxldmVscyA9IGMoXCJGQUxTRVwiLCBcIlRSVUVcIikpKVxuICBcbiAgbG9oX2RmW2dlbmUsIGMoXCJuX25vX3BscF9ub19hbHRcIiwgXCJuX3BscF9ub19hbHRcIiwgXCJuX25vX3BscF9hbHRcIiwgXCJuX3BscF9hbHRcIildID0gdW5saXN0KGxvaF9tYXQpXG4gIGxvaF9kZltnZW5lLCBjKFwib2Rkc1JhdGlvXCIsIFwicFwiLCBcIkNJX2xvd2VyXCIsIFwiQ0lfdXBwZXJcIildIDwtIHVubGlzdChmaXNoZXIudGVzdChsb2hfbWF0KVtjKFwiZXN0aW1hdGVcIiwgXCJwLnZhbHVlXCIsIFwiY29uZi5pbnRcIildKVxuXG4gICMgRmlsdGVyIGV4cHIgZGF0YSBmb3Igc2FtcGxlcyB3aXRoIGV4cHJlc3Npb24gbG9zcyBpbiBnZW5lXG4gIGdlbmVfZXhwcl9sb3NzIDwtIHpfc2NvcmVkX2V4cHJlc3Npb24gJT4lXG4gICAgZHBseXI6OmZpbHRlcihIdWdvX1N5bWJvbCA9PSBnZW5lICYgZXhwcl96c2NvcmUgPCAtMilcbiAgXG4gICMgY3JlYXRlIFAvTFAgYW5kIGV4cHIgbG9zcyBjb250aW5nZW5jeSB0YWJsZSBhbmQgc3RvcmUgc3RhdHMgaW4gZXhwcl9sb3NzX2RmXG4gIGV4cHJfbG9zc19tYXQgPC0gdGFibGUoaXNfcGxwID0gaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxbaGlzdF9jb2hvcnQkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwgJWluJSB6X3Njb3JlZF9leHByZXNzaW9uJEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsXSAlaW4lIHBscF9nZW5lX2JzLFxuICAgICAgICAgICAgICAgICAgIGlzX2FsdCA9IGZhY3RvcihoaXN0X2NvaG9ydCRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbFtoaXN0X2NvaG9ydCRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCAlaW4lIHpfc2NvcmVkX2V4cHJlc3Npb24kS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxdICVpbiUgZ2VuZV9leHByX2xvc3MkS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWwsIGxldmVscyA9IGMoXCJGQUxTRVwiLCBcIlRSVUVcIikpKVxuICBcbiAgZXhwcl9sb3NzX2RmW2dlbmUsIGMoXCJuX25vX3BscF9ub19hbHRcIiwgXCJuX3BscF9ub19hbHRcIiwgXCJuX25vX3BscF9hbHRcIiwgXCJuX3BscF9hbHRcIildID0gdW5saXN0KGV4cHJfbG9zc19tYXQpXG4gIGV4cHJfbG9zc19kZltnZW5lLCBjKFwib2Rkc1JhdGlvXCIsIFwicFwiLCBcIkNJX2xvd2VyXCIsIFwiQ0lfdXBwZXJcIildIDwtIHVubGlzdChmaXNoZXIudGVzdChleHByX2xvc3NfbWF0KVtjKFwiZXN0aW1hdGVcIiwgXCJwLnZhbHVlXCIsIFwiY29uZi5pbnRcIildKVxuXG4gICMgUmVwZWF0IGZvciBleHByZXNzaW9uIGdhaW5zXG4gIGdlbmVfZXhwcl9nYWluIDwtIHpfc2NvcmVkX2V4cHJlc3Npb24gJT4lXG4gICAgZHBseXI6OmZpbHRlcihIdWdvX1N5bWJvbCA9PSBnZW5lICYgZXhwcl96c2NvcmUgPiAyKVxuICBcbiAgZXhwcl9nYWluX21hdCA8LSB0YWJsZShpc19wbHAgPSBoaXN0X2NvaG9ydCRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbFtoaXN0X2NvaG9ydCRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCAlaW4lIHpfc2NvcmVkX2V4cHJlc3Npb24kS2lkc19GaXJzdF9CaW9zcGVjaW1lbl9JRF9ub3JtYWxdICVpbiUgcGxwX2dlbmVfYnMsXG4gICAgICAgICAgICAgICAgICAgaXNfYWx0ID0gZmFjdG9yKGhpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsW2hpc3RfY29ob3J0JEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfbm9ybWFsICVpbiUgel9zY29yZWRfZXhwcmVzc2lvbiRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbF0gJWluJSBnZW5lX2V4cHJfZ2FpbiRLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCwgbGV2ZWxzID0gYyhcIkZBTFNFXCIsIFwiVFJVRVwiKSkpXG4gIFxuICBleHByX2dhaW5fZGZbZ2VuZSwgYyhcIm5fbm9fcGxwX25vX2FsdFwiLCBcIm5fcGxwX25vX2FsdFwiLCBcIm5fbm9fcGxwX2FsdFwiLCBcIm5fcGxwX2FsdFwiKV0gPSB1bmxpc3QoZXhwcl9nYWluX21hdClcbiAgZXhwcl9nYWluX2RmW2dlbmUsIGMoXCJvZGRzUmF0aW9cIiwgXCJwXCIsIFwiQ0lfbG93ZXJcIiwgXCJDSV91cHBlclwiKV0gPC0gdW5saXN0KGZpc2hlci50ZXN0KGV4cHJfZ2Fpbl9tYXQpW2MoXCJlc3RpbWF0ZVwiLCBcInAudmFsdWVcIiwgXCJjb25mLmludFwiKV0pXG5cblxufSAgIn0= -->

```r
for (gene in genes){
  
  # identify samples with P-LP variant in gene
  plp_gene_bs <- plp %>%
    dplyr::filter(Hugo_Symbol == gene) %>%
    pull(Kids_First_Biospecimen_ID) 

  # Filter MAF for only variants in gene
  gene_maf <- maf %>%
      dplyr::filter(Hugo_Symbol == gene & ONCOGENIC %in% c("Oncogenic", "Likely Oncogenic")) 
  
  # create contingency table of plp T/F and somatic alteration T/F
  snv_mat <- table(is_plp = hist_cohort$Kids_First_Biospecimen_ID_normal %in% plp_gene_bs, 
                   is_alt = factor(hist_cohort$Kids_First_Biospecimen_ID_normal %in% gene_maf$Matched_Norm_Sample_Barcode, levels = c("FALSE", "TRUE")))
  
  # add stats to df 
  snv_df[gene, c("n_no_plp_no_alt", "n_plp_no_alt", "n_no_plp_alt", "n_plp_alt")] = unlist(snv_mat)
  snv_df[gene, c("oddsRatio", "p", "CI_lower", "CI_upper")] <- unlist(fisher.test(snv_mat)[c("estimate", "p.value", "conf.int")])
  
  # filter cnv df for samples with loss in gene
  gene_cnv_loss <- cnv %>%
    dplyr::filter(gene_symbol == gene & status == "loss")
  
  # create P/LP and somatic alteration contingency table and store stats in df
  cnv_loss_mat <- table(is_plp = hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% cnv$Kids_First_Biospecimen_ID_normal] %in% plp_gene_bs,
                        is_alt = factor(hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% cnv$Kids_First_Biospecimen_ID_normal] %in% gene_cnv_loss$Kids_First_Biospecimen_ID_normal, levels = c("FALSE", "TRUE")))
  
  cnv_loss_df[gene, c("n_no_plp_no_alt", "n_plp_no_alt", "n_no_plp_alt", "n_plp_alt")] = unlist(cnv_loss_mat)
  cnv_loss_df[gene, c("oddsRatio", "p", "CI_lower", "CI_upper")] <- unlist(fisher.test(cnv_loss_mat)[c("estimate", "p.value", "conf.int")])
  
  # Repeat for copy number gains
  gene_cnv_gain <- cnv %>%
    dplyr::filter(gene_symbol == gene & status == "gain")
  
  cnv_gain_mat <- table(is_plp = hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% cnv$Kids_First_Biospecimen_ID_normal] %in% plp_gene_bs,
                        is_alt = factor(hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% cnv$Kids_First_Biospecimen_ID_normal] %in% gene_cnv_gain$Kids_First_Biospecimen_ID_normal, levels = c("FALSE", "TRUE")))
  
  cnv_gain_df[gene, c("n_no_plp_no_alt", "n_plp_no_alt", "n_no_plp_alt", "n_plp_alt")] = unlist(cnv_gain_mat)
  cnv_gain_df[gene, c("oddsRatio", "p", "CI_lower", "CI_upper")] <- unlist(fisher.test(cnv_gain_mat)[c("estimate", "p.value", "conf.int")])
  
  # Filter LOH df for samples exhibiting LOH in gene
  gene_loh <- loh %>%
    dplyr::filter(Hugo_Symbol == gene & LOH_score > 0.25)
  
  # Create contingency table and store stats for LOH enrichment 
  loh_mat <- table(is_plp = hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% loh$Kids_First_Biospecimen_ID_normal] %in% plp_gene_bs,
                        is_alt = factor(hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% loh$Kids_First_Biospecimen_ID_normal] %in% gene_loh$Kids_First_Biospecimen_ID_normal, levels = c("FALSE", "TRUE")))
  
  loh_df[gene, c("n_no_plp_no_alt", "n_plp_no_alt", "n_no_plp_alt", "n_plp_alt")] = unlist(loh_mat)
  loh_df[gene, c("oddsRatio", "p", "CI_lower", "CI_upper")] <- unlist(fisher.test(loh_mat)[c("estimate", "p.value", "conf.int")])

  # Filter expr data for samples with expression loss in gene
  gene_expr_loss <- z_scored_expression %>%
    dplyr::filter(Hugo_Symbol == gene & expr_zscore < -2)
  
  # create P/LP and expr loss contingency table and store stats in expr_loss_df
  expr_loss_mat <- table(is_plp = hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% z_scored_expression$Kids_First_Biospecimen_ID_normal] %in% plp_gene_bs,
                   is_alt = factor(hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% z_scored_expression$Kids_First_Biospecimen_ID_normal] %in% gene_expr_loss$Kids_First_Biospecimen_ID_normal, levels = c("FALSE", "TRUE")))
  
  expr_loss_df[gene, c("n_no_plp_no_alt", "n_plp_no_alt", "n_no_plp_alt", "n_plp_alt")] = unlist(expr_loss_mat)
  expr_loss_df[gene, c("oddsRatio", "p", "CI_lower", "CI_upper")] <- unlist(fisher.test(expr_loss_mat)[c("estimate", "p.value", "conf.int")])

  # Repeat for expression gains
  gene_expr_gain <- z_scored_expression %>%
    dplyr::filter(Hugo_Symbol == gene & expr_zscore > 2)
  
  expr_gain_mat <- table(is_plp = hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% z_scored_expression$Kids_First_Biospecimen_ID_normal] %in% plp_gene_bs,
                   is_alt = factor(hist_cohort$Kids_First_Biospecimen_ID_normal[hist_cohort$Kids_First_Biospecimen_ID_normal %in% z_scored_expression$Kids_First_Biospecimen_ID_normal] %in% gene_expr_gain$Kids_First_Biospecimen_ID_normal, levels = c("FALSE", "TRUE")))
  
  expr_gain_df[gene, c("n_no_plp_no_alt", "n_plp_no_alt", "n_no_plp_alt", "n_plp_alt")] = unlist(expr_gain_mat)
  expr_gain_df[gene, c("oddsRatio", "p", "CI_lower", "CI_upper")] <- unlist(fisher.test(expr_gain_mat)[c("estimate", "p.value", "conf.int")])


}  
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Perform multiple test correction for each alteration type

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbInNudl9kZiRhZGp1c3RlZF9wIDwtIHAuYWRqdXN0KHNudl9kZiRwLCBtZXRob2QgPSBcIkJIXCIpIiwiY252X2xvc3NfZGYkYWRqdXN0ZWRfcCA8LSBwLmFkanVzdChjbnZfbG9zc19kZiRwLCBtZXRob2QgPSBcIkJIXCIpIiwiY252X2dhaW5fZGYkYWRqdXN0ZWRfcCA8LSBwLmFkanVzdChjbnZfZ2Fpbl9kZiRwLCBtZXRob2QgPSBcIkJIXCIpIiwibG9oX2RmJGFkanVzdGVkX3AgPC0gcC5hZGp1c3QobG9oX2RmJHAsIG1ldGhvZCA9IFwiQkhcIikiLCJleHByX2xvc3NfZGYkYWRqdXN0ZWRfcCA8LSBwLmFkanVzdChleHByX2xvc3NfZGYkcCwgbWV0aG9kID0gXCJCSFwiKSIsImV4cHJfZ2Fpbl9kZiRhZGp1c3RlZF9wIDwtIHAuYWRqdXN0KGV4cHJfZ2Fpbl9kZiRwLCBtZXRob2QgPSBcIkJIXCIpIl19 -->

```r
snv_df$adjusted_p <- p.adjust(snv_df$p, method = "BH")
cnv_loss_df$adjusted_p <- p.adjust(cnv_loss_df$p, method = "BH")
cnv_gain_df$adjusted_p <- p.adjust(cnv_gain_df$p, method = "BH")
loh_df$adjusted_p <- p.adjust(loh_df$p, method = "BH")
expr_loss_df$adjusted_p <- p.adjust(expr_loss_df$p, method = "BH")
expr_gain_df$adjusted_p <- p.adjust(expr_gain_df$p, method = "BH")
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Create enrichment, p-value, logical significance, and label matrices for plotting

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbImFsdF9kZiA8LSBzbnZfZGYgJT4lXG4gIGJpbmRfcm93cyhjbnZfbG9zc19kZiwgY252X2dhaW5fZGYsIFxuICAgICAgICAgICAgbG9oX2RmLCBleHByX2xvc3NfZGYsIGV4cHJfZ2Fpbl9kZikgJT4lXG4gIGRwbHlyOjptdXRhdGUoSHVnb19TeW1ib2wgPSByZXAoZ2VuZXMsIDYpKSIsIiIsImFsdGVyYXRpb25zID0gYyhcIk9uY29nZW5pYyBTTlZcIiwgXCJDTiBsb3NzXCIsIFwiQ04gZ2FpblwiLFxuICAgICAgICAgICAgICAgIFwiTE9IID4gMC4yNVwiLCBcIkV4cHIgVFBNIDwgLTJcIiwgXCJFeHByIFRQTSA+IDJcIikiLCIiLCJlbnJfbWF0IDwtIG1hdHJpeChhbHRfZGYkb2Rkc1JhdGlvLCBcbiAgICAgICAgICAgICAgICAgIGxlbmd0aChnZW5lcyksIGxlbmd0aChhbHRlcmF0aW9ucyksXG4gICAgICAgICAgICAgICAgICBkaW1uYW1lcyA9IGxpc3QoZ2VuZXMsIGFsdGVyYXRpb25zKSkiLCIiLCJwX21hdCA8LSBtYXRyaXgoYWx0X2RmJGFkanVzdGVkX3AsIFxuICAgICAgICAgICAgICAgICAgbGVuZ3RoKGdlbmVzKSwgbGVuZ3RoKGFsdGVyYXRpb25zKSxcbiAgICAgICAgICAgICAgICAgIGRpbW5hbWVzID0gbGlzdChnZW5lcywgYWx0ZXJhdGlvbnMpKSIsIiIsInNpZ19tYXQgPC0gaWZlbHNlKHBfbWF0IDwgMC4wNSwgXCIqXCIsIFwiXCIpIiwiIiwiZmlsbF9tYXQgPC0gbWF0cml4KGdsdWU6OmdsdWUoXCJ7cm91bmQoZW5yX21hdCwgMSl9IHtzaWdfbWF0fVwiKSwgXG4gICAgICAgICAgICAgICAgICAgbGVuZ3RoKGdlbmVzKSwgbGVuZ3RoKGFsdGVyYXRpb25zKSkiXX0= -->

```r
alt_df <- snv_df %>%
  bind_rows(cnv_loss_df, cnv_gain_df, 
            loh_df, expr_loss_df, expr_gain_df) %>%
  dplyr::mutate(Hugo_Symbol = rep(genes, 6))

alterations = c("Oncogenic SNV", "CN loss", "CN gain",
                "LOH > 0.25", "Expr TPM < -2", "Expr TPM > 2")

enr_mat <- matrix(alt_df$oddsRatio, 
                  length(genes), length(alterations),
                  dimnames = list(genes, alterations))

p_mat <- matrix(alt_df$adjusted_p, 
                  length(genes), length(alterations),
                  dimnames = list(genes, alterations))

sig_mat <- ifelse(p_mat < 0.05, "*", "")

fill_mat <- matrix(glue::glue("{round(enr_mat, 1)} {sig_mat}"), 
                   length(genes), length(alterations))
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Plot somatic alteration enrichment heatmap

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbIiMgUGxvdCBlbnJpY2htZW50IGhlYXRtYXAiLCJjb2xfZnVuID0gY29sb3JSYW1wMihjKDAsIDYpLCBjKFwid2hpdGVcIiwgXCJvcmFuZ2VyZWRcIikpIiwiIiwicGRmKGZpbGUucGF0aChwbG90X2RpciwgXCJzb21hdGljLWFsdGVyYXRpb24tZW5yLWhlYXRtYXAucGRmXCIpLFxuICAgICAgd2lkdGggPSA1LCBoZWlnaHQgPSA2KSIsIiIsImVucl9oZWF0bWFwIDwtIEhlYXRtYXAoZW5yX21hdCxcbiAgICAgICAgICAgICAgICAgICAgICAgbmFtZSA9IFwiT2RkcyByYXRpb1wiLFxuICAgICAgICAgICAgICAgICAgICAgICBjbHVzdGVyX3Jvd3MgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICBjbHVzdGVyX2NvbHVtbnMgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICByZWN0X2dwID0gZ3Bhcihjb2wgPSBcImJsYWNrXCIsIGx3ZCA9IDIpLFxuICAgICAgICAgICAgICAgICAgICAgICBjb2wgPSBjb2xfZnVuLFxuICAgICAgICAgICAgICAgICAgICAgICBjZWxsX2Z1biA9IGZ1bmN0aW9uKGosIGksIHgsIHksIHdpZHRoLCBoZWlnaHQsIGZpbGwpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgZ3JpZC50ZXh0KHNwcmludGYoXCIlc1wiLCBmaWxsX21hdFtpLCBqXSksIHgsIHksIGdwID0gZ3Bhcihmb250c2l6ZSA9IDEwKSlcbiAgICAgICAgICAgICAgIH0pIiwiIiwicHJpbnQoZW5yX2hlYXRtYXApIiwiIiwiZGV2Lm9mZigpIl19 -->

```r
# Plot enrichment heatmap
col_fun = colorRamp2(c(0, 6), c("white", "orangered"))

pdf(file.path(plot_dir, "somatic-alteration-enr-heatmap.pdf"),
      width = 5, height = 6)

enr_heatmap <- Heatmap(enr_mat,
                       name = "Odds ratio",
                       cluster_rows = F,
                       cluster_columns = F,
                       rect_gp = gpar(col = "black", lwd = 2),
                       col = col_fun,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(sprintf("%s", fill_mat[i, j]), x, y, gp = gpar(fontsize = 10))
               })

print(enr_heatmap)

dev.off()
```



<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoicG5nIFxuICAyIFxuIn0= -->

```
png 
  2 
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Create data frame of significant enrichments to plot

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbInNpZ19hbHRfZGYgPC0gYWx0X2RmICU+JVxuICBkcGx5cjo6ZmlsdGVyKGFkanVzdGVkX3AgPCAwLjA1KSAlPiVcbiAgZHBseXI6OnNlbGVjdChIdWdvX1N5bWJvbCwgLW5fbm9fcGxwX2FsdCwgLW5fbm9fcGxwX25vX2FsdCwgZXZlcnl0aGluZygpKSAlPiVcbiAgZHBseXI6OnJlbmFtZShuX2FsdCA9IG5fcGxwX2FsdCxcbiAgICAgICAgICAgICAgICBuX25vX2FsdCA9IG5fcGxwX25vX2FsdCkgJT4lXG4gIGRwbHlyOjptdXRhdGUoZ3JvdXAgPSBcIkdlcm1saW5lIFAtTFBcIikiLCIiLCIjIGNyZWF0ZSByZWZlcmVuY2UgZGF0YSBmcmFtZSB3aXRoIGVucmljaG1lbnQgPSAxIHRvIHBsb3QgYWdhaW5zdCBlbnJpY2htZW50IHNjb3JlcyBpbiBwYXRpZW50cyB3LyBQLUxQIHZhcmlhbnRzIiwiZGZfbm9fcGxwIDwtIHNpZ19hbHRfZGYgJT4lXG4gIGRwbHlyOjpzZWxlY3QoSHVnb19TeW1ib2wsIGFsdGVyYXRpb24sIFxuICAgICAgICAgICAgICAgIG5fbm9fcGxwX2FsdCwgbl9ub19wbHBfbm9fYWx0KSAlPiVcbiAgZHBseXI6OnJlbmFtZShuX2FsdCA9IG5fbm9fcGxwX2FsdCxcbiAgICAgICAgICAgICAgICBuX25vX2FsdCA9IG5fbm9fcGxwX25vX2FsdCkgJT4lXG4gIGRwbHlyOjptdXRhdGUob2Rkc1JhdGlvID0gMSxcbiAgICAgICAgICAgICAgICBwID0gTkEsIFxuICAgICAgICAgICAgICAgIENJX2xvd2VyID0gTkEsIENJX3VwcGVyID0gTkEsIFxuICAgICAgICAgICAgICAgIGFkanVzdGVkX3AgPSBOQSxcbiAgICAgICAgICAgICAgICBncm91cCA9IFwiTm8gZ2VybWxpbmUgUC1MUFwiKSIsIiIsIiMgbWVyZ2UgZGF0YSBmcmFtZXMsIGFuZCBjcmVhdGUgYWRkaXRpb25hbCBjb2x1bW5zIGZvciBwbG90dGluZyIsImRmX2FsbCA8LSBzaWdfYWx0X2RmICU+JVxuICBiaW5kX3Jvd3MoZGZfbm9fcGxwKSAlPiVcbiAgZHBseXI6Om11dGF0ZShwZXJjID0gbl9hbHQvKG5fYWx0ICsgbl9ub19hbHQpICogMTAwLFxuICAgICAgICAgICAgICAgIGZyYWN0aW9uID0gZ2x1ZTo6Z2x1ZShcIntuX2FsdH0ve25fYWx0ICsgbl9ub19hbHR9XCIpLFxuICAgICAgICAgICAgICAgIHBfdGV4dCA9IGNhc2Vfd2hlbihcbiAgICAgICAgICAgICAgICAgIGdyb3VwID09IFwiR2VybWxpbmUgUC1MUFwiIH4gZ2x1ZTo6Z2x1ZShcInA9e3NpZ25pZihhZGp1c3RlZF9wLCBkaWdpdHMgPSAyKX1cIiksXG4gICAgICAgICAgICAgICAgICBUUlVFIH4gTkFfY2hhcmFjdGVyX1xuICAgICAgICAgICAgICAgICkpICU+JVxuICBkcGx5cjo6bXV0YXRlKGdlbmVfYWx0ID0gZ2x1ZTo6Z2x1ZShcIntIdWdvX1N5bWJvbH06IHthbHRlcmF0aW9ufVwiKSkiXX0= -->

```r
sig_alt_df <- alt_df %>%
  dplyr::filter(adjusted_p < 0.05) %>%
  dplyr::select(Hugo_Symbol, -n_no_plp_alt, -n_no_plp_no_alt, everything()) %>%
  dplyr::rename(n_alt = n_plp_alt,
                n_no_alt = n_plp_no_alt) %>%
  dplyr::mutate(group = "Germline P-LP")

# create reference data frame with enrichment = 1 to plot against enrichment scores in patients w/ P-LP variants
df_no_plp <- sig_alt_df %>%
  dplyr::select(Hugo_Symbol, alteration, 
                n_no_plp_alt, n_no_plp_no_alt) %>%
  dplyr::rename(n_alt = n_no_plp_alt,
                n_no_alt = n_no_plp_no_alt) %>%
  dplyr::mutate(oddsRatio = 1,
                p = NA, 
                CI_lower = NA, CI_upper = NA, 
                adjusted_p = NA,
                group = "No germline P-LP")

# merge data frames, and create additional columns for plotting
df_all <- sig_alt_df %>%
  bind_rows(df_no_plp) %>%
  dplyr::mutate(perc = n_alt/(n_alt + n_no_alt) * 100,
                fraction = glue::glue("{n_alt}/{n_alt + n_no_alt}"),
                p_text = case_when(
                  group == "Germline P-LP" ~ glue::glue("p={signif(adjusted_p, digits = 2)}"),
                  TRUE ~ NA_character_
                )) %>%
  dplyr::mutate(gene_alt = glue::glue("{Hugo_Symbol}: {alteration}"))
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Plot enrichment of somatic alterations in patients with P-LP variants relative to those without

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbIiMgQ3JlYXRlIENQRyBPZGRzIFJhdGlvIHBsb3QgIiwiZW5yX3Bsb3QgPC0gZGZfYWxsICU+JSBcbiAgZ2dwbG90KGFlcyh4ID0gZ3JvdXAsIHkgPSBsb2cxMChvZGRzUmF0aW8pLCBsYWJlbCA9IHBfdGV4dCkpICtcbiAgZ2VvbV9wb2ludChzaXplID0gMywgY29sb3IgPSBcIiMwMEEwODdGRlwiLFxuICAgICAgICAgICAgIHNob3cubGVnZW5kID0gRkFMU0UpICsgXG4gIGdlb21fZXJyb3JiYXIoYWVzKHltaW4gPSBsb2cxMChDSV9sb3dlciksIHltYXggPSBsb2cxMChDSV91cHBlcikpLCB3aWR0aCA9IDAuMiwgXG4gICAgICAgICAgICAgICAgc2hvdy5sZWdlbmQgPSBGQUxTRSwgY29sb3IgPSBcIiMwMEEwODdGRlwiKSArXG4gIGxhYnMoeSA9IFwibG9nMTAtT1IgKDk1JSBDSSlcIiwgeCA9IE5VTEwpICsgXG4gIGdlb21fdGV4dCh5ID0gMS42LCB4ID0gMiwgaGp1c3QgPSAwLCBzaXplID0gNCwgZm9udGZhY2UgPSAyKSArXG4gIGNvb3JkX2ZsaXAoKSArXG4gIGdlb21faGxpbmUoeWludGVyY2VwdCA9IGxvZzEwKDEpLCBsaW5ldHlwZSA9IFwiZGFzaGVkXCIpICsgXG4gIGZhY2V0X3dyYXAofmdlbmVfYWx0LCBucm93ID0gNywgc2NhbGUgPSBcImZpeGVkXCIpICtcbiAgZXhwYW5kX2xpbWl0cyh5PTApICtcbiAgdGhlbWVfUHVibGljYXRpb24oKSArXG4gIHRoZW1lKHBsb3QubWFyZ2luID0gdW5pdChjKDIsMC41LDEsMC4yNSksIFwibGluZXNcIikpIiwiIiwiIyBDcmVhdGUgJSBwYXRpZW50cyB3aXRoIENQRyBQTFAgcGxvdCBzZXBhcmF0ZWQgYnkgc291cmNlIGNhbGwsIGFuZCBpbmNsdWRlIGZyYWN0aW9ucyBhcyB0ZXh0ICIsInBlcmNfcGxvdCA8LSBkZl9hbGwgJT4lXG4gIGdncGxvdChhZXMoeCA9IHBlcmMsIHkgPSBncm91cCwgbGFiZWwgPSBmcmFjdGlvbiwgZmlsbCA9IGdyb3VwKSkgK1xuICBnZW9tX2JhcihzdGF0ID0gXCJpZGVudGl0eVwiLCBjb2xvciA9IFwiYmxhY2tcIixcbiAgICAgICAgICAgc2hvdy5sZWdlbmQgPSBGQUxTRSkgKyBcbiAgZ2VvbV90ZXh0KHggPSA4NSwgaGp1c3QgPSAwLCBzaXplID0gNCwgZm9udGZhY2UgPSAyKSArXG4gIGxhYnMoeCA9IFwiJSBwYXRpZW50c1wiLCB5ID0gTlVMTCwgZmlsbCA9IE5VTEwpICsgXG4gIGd1aWRlcyhmaWxsID0gZ3VpZGVfbGVnZW5kKG5yb3cgPSAxKSkgK1xuICBmYWNldF93cmFwKH5nZW5lX2FsdCwgbnJvdyA9IDcsIHNjYWxlID0gXCJmaXhlZFwiKSArXG4gIHNjYWxlX3lfZGlzY3JldGUobGFiZWxzID0gYyhcIkdlcm1saW5lIFAtTFBcIiA9IE5VTEwsXG4gICAgICAgICAgICAgICAgICAgXCJObyBnZXJtbGluZSBQLUxQXCIgPSBOVUxMKSkgK1xuICBleHBhbmRfbGltaXRzKHg9MikgK1xuICBjb29yZF9jYXJ0ZXNpYW4oY2xpcCA9ICdvZmYnKSArXG4gIHRoZW1lX1B1YmxpY2F0aW9uKCkgK1xuICB0aGVtZShwbG90Lm1hcmdpbiA9IHVuaXQoYygyLDYsMSwxKSwgXCJsaW5lc1wiKSxcbiAgICAgICAgbGVnZW5kLnBvc2l0aW9uID0gYygwLjUsIDEuMDcpKSIsIiIsIiIsInBkZihmaWxlLnBhdGgocGxvdF9kaXIsIFwic2lnLXNvbWF0aWMtYWx0ZXJhdGlvbi1lbnIucGRmXCIpLFxuICAgICB3aWR0aCA9IDcsIGhlaWdodCA9IDgpIiwiIiwiZ2dhcnJhbmdlKGVucl9wbG90LCBwZXJjX3Bsb3QsXG4gICAgICAgICAgbnJvdyA9IDEsIHdpZHRocyA9IGMoMiwyKSkiXX0= -->

```r
# Create CPG Odds Ratio plot 
enr_plot <- df_all %>% 
  ggplot(aes(x = group, y = log10(oddsRatio), label = p_text)) +
  geom_point(size = 3, color = "#00A087FF",
             show.legend = FALSE) + 
  geom_errorbar(aes(ymin = log10(CI_lower), ymax = log10(CI_upper)), width = 0.2, 
                show.legend = FALSE, color = "#00A087FF") +
  labs(y = "log10-OR (95% CI)", x = NULL) + 
  geom_text(y = 1.6, x = 2, hjust = 0, size = 4, fontface = 2) +
  coord_flip() +
  geom_hline(yintercept = log10(1), linetype = "dashed") + 
  facet_wrap(~gene_alt, nrow = 7, scale = "fixed") +
  expand_limits(y=0) +
  theme_Publication() +
  theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"))

# Create % patients with CPG PLP plot separated by source call, and include fractions as text 
perc_plot <- df_all %>%
  ggplot(aes(x = perc, y = group, label = fraction, fill = group)) +
  geom_bar(stat = "identity", color = "black",
           show.legend = FALSE) + 
  geom_text(x = 85, hjust = 0, size = 4, fontface = 2) +
  labs(x = "% patients", y = NULL, fill = NULL) + 
  guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(~gene_alt, nrow = 7, scale = "fixed") +
  scale_y_discrete(labels = c("Germline P-LP" = NULL,
                   "No germline P-LP" = NULL)) +
  expand_limits(x=2) +
  coord_cartesian(clip = 'off') +
  theme_Publication() +
  theme(plot.margin = unit(c(2,6,1,1), "lines"),
        legend.position = c(0.5, 1.07))


pdf(file.path(plot_dir, "sig-somatic-alteration-enr.pdf"),
     width = 7, height = 8)

ggarrange(enr_plot, perc_plot,
          nrow = 1, widths = c(2,2))
```



<!-- rnb-source-end -->

<!-- rnb-warning-begin eyJkYXRhIjoiV2FybmluZzogUmVtb3ZlZCA1IHJvd3MgY29udGFpbmluZyBtaXNzaW5nIHZhbHVlcyAoYGdlb21fdGV4dCgpYCkuXG4ifQ== -->

```
Warning: Removed 5 rows containing missing values (`geom_text()`).
```



<!-- rnb-warning-end -->

<!-- rnb-source-begin eyJkYXRhIjoiZGV2Lm9mZigpIn0= -->

```r
dev.off()
```



<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoicG5nIFxuICAyIFxuIn0= -->

```
png 
  2 
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Write `alt_df` to output

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYWx0X2RmICU+JVxuICBkcGx5cjo6cmVsb2NhdGUoSHVnb19TeW1ib2wpICU+JVxuICB3cml0ZV90c3YoZmlsZS5wYXRoKHJlc3VsdHNfZGlyLCBcInNvbWF0aWMtYWx0ZXJhdGlvbi1lbnJpY2htZW50LnRzdlwiKSkifQ== -->

```r
alt_df %>%
  dplyr::relocate(Hugo_Symbol) %>%
  write_tsv(file.path(results_dir, "somatic-alteration-enrichment.tsv"))
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->

