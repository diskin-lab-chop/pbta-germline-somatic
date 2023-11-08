
<!-- rnb-text-begin -->

---
title: "Annotate MAF"
output: 
  html_notebook:
  toc: TRUE
toc_float: TRUE
author: Jo Lynne Rokita
params:
  maf_in:
    label: "MAF input file"
    value: data/snv-consensus-plus-hotspots.maf.tsv.gz
    input: file
  maf_goi:
    label: "goi MAF output file"
    value: results/snv-consensus-plus-hotspots-goi.maf.tsv
    input: file
---


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbIiMgbG9hZCBsaWJyYXJpZXMiLCJsaWJyYXJ5KHRpZHl2ZXJzZSkiXX0= -->

```r
# load libraries
library(tidyverse)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoi4pSA4pSAIEF0dGFjaGluZyBwYWNrYWdlcyDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIDilIAgdGlkeXZlcnNlIDEuMy4yIOKUgOKUgFxu4pyUIGdncGxvdDIgMy40LjAgICAgIOKclCBwdXJyciAgIDEuMC4xXG7inJQgdGliYmxlICAzLjEuOCAgICAg4pyUIGRwbHlyICAgMS4xLjBcbuKclCB0aWR5ciAgIDEuMy4wICAgICDinJQgc3RyaW5nciAxLjUuMFxu4pyUIHJlYWRyICAgMi4xLjMgICAgIOKclCBmb3JjYXRzIDEuMC4wXG7ilIDilIAgQ29uZmxpY3RzIOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgOKUgCB0aWR5dmVyc2VfY29uZmxpY3RzKCkg4pSA4pSAXG7inJYgZHBseXI6OmZpbHRlcigpIG1hc2tzIHN0YXRzOjpmaWx0ZXIoKVxu4pyWIGRwbHlyOjpsYWcoKSAgICBtYXNrcyBzdGF0czo6bGFnKClcbiJ9 -->

```
── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.2 ──
✔ ggplot2 3.4.0     ✔ purrr   1.0.1
✔ tibble  3.1.8     ✔ dplyr   1.1.0
✔ tidyr   1.3.0     ✔ stringr 1.5.0
✔ readr   2.1.3     ✔ forcats 1.0.0
── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjpbIiMgZGVmaW5lIGRpcmVjdG9yaWVzICIsInJvb3RfZGlyIDwtIHJwcm9qcm9vdDo6ZmluZF9yb290KHJwcm9qcm9vdDo6aGFzX2RpcihcIi5naXRcIikpIiwiYW5hbHlzaXNfZGlyIDwtIGZpbGUucGF0aChyb290X2RpciwgXCJhbmFseXNlc1wiLCBcIm9uY29rYi1hbm5vdGF0aW9uXCIpIiwiaW5wdXRfZGlyIDwtIGZpbGUucGF0aChhbmFseXNpc19kaXIsIFwiaW5wdXRcIikiLCJyZXN1bHRzX2RpciA8LSBmaWxlLnBhdGgoYW5hbHlzaXNfZGlyLCBcInJlc3VsdHNcIikiLCJkYXRhX2RpciA8LSBmaWxlLnBhdGgocm9vdF9kaXIsIFwiZGF0YVwiKSIsIiIsIiMgSWYgdGhlIGlucHV0IGFuZCByZXN1bHRzIGRpcmVjdG9yaWVzIGRvIG5vdCBleGlzdCwgY3JlYXRlIGl0IiwiaWYgKCFkaXIuZXhpc3RzKHJlc3VsdHNfZGlyKSkge1xuICBkaXIuY3JlYXRlKHJlc3VsdHNfZGlyLCByZWN1cnNpdmUgPSBUUlVFKVxufSIsIiIsImlmICghZGlyLmV4aXN0cyhpbnB1dF9kaXIpKSB7XG4gIGRpci5jcmVhdGUoaW5wdXRfZGlyLCByZWN1cnNpdmUgPSBUUlVFKVxufSIsIiIsIiMgYWRkIGhpc3RvbG9naWVzIGZpbGUgc28gd2UgY2FuIHN1YnNldCBtYWYiLCJoaXN0b2xvZ2llc19maWxlIDwtIGZpbGUucGF0aChyb290X2RpciwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICBcImFuYWx5c2VzXCIsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjb2xsYXBzZS10dW1vci1oaXN0b2xvZ2llc1wiLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgIFwicmVzdWx0c1wiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJnZXJtbGluZS1wcmltYXJ5LXBsdXMtdHVtb3ItaGlzdG9sb2dpZXMtcGxvdC1ncm91cHMtY2xpbi1tZXRhLnRzdlwiKSIsIiIsIiMgZ29pIGZpbGUiLCJjcGdfZmlsZSA8LSBmaWxlLnBhdGgoYW5hbHlzaXNfZGlyLCBcImlucHV0XCIsIFwiY3BnLnR4dFwiKSJdfQ== -->

```r
# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "oncokb-annotation")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
data_dir <- file.path(root_dir, "data")

# If the input and results directories do not exist, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE)
}

# add histologies file so we can subset maf
histologies_file <- file.path(root_dir, 
                           "analyses", 
                           "collapse-tumor-histologies", 
                           "results",
                           "germline-primary-plus-tumor-histologies-plot-groups-clin-meta.tsv")

# goi file
cpg_file <- file.path(analysis_dir, "input", "cpg.txt")
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoidHVtb3JfaWRzIDwtIHJlYWRfdHN2KGhpc3RvbG9naWVzX2ZpbGUsIGd1ZXNzX21heCA9IDMwMDApICU+JVxuICBwdWxsKEtpZHNfRmlyc3RfQmlvc3BlY2ltZW5fSURfdHVtb3IpIn0= -->

```r
tumor_ids <- read_tsv(histologies_file, guess_max = 3000) %>%
  pull(Kids_First_Biospecimen_ID_tumor)
```



<!-- rnb-source-end -->

<!-- rnb-message-begin eyJkYXRhIjoiUm93czogODM4IENvbHVtbnM6IDM2XG7ilIDilIAgQ29sdW1uIHNwZWNpZmljYXRpb24g4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSAXG5EZWxpbWl0ZXI6IFwiXFx0XCJcbmNociAoMjYpOiBLaWRzX0ZpcnN0X1BhcnRpY2lwYW50X0lELCBLaWRzX0ZpcnN0X0Jpb3NwZWNpbWVuX0lEX25vcm1hbCwgc2FtcGwuLi5cbmRibCAoMTApOiBicm9hZF9oaXN0b2xvZ3lfb3JkZXIsIGFnZV9hdF9kaWFnbm9zaXNfZGF5cywgYWdlX2xhc3RfdXBkYXRlX2RheXMuLi5cblxu4oS5IFVzZSBgc3BlYygpYCB0byByZXRyaWV2ZSB0aGUgZnVsbCBjb2x1bW4gc3BlY2lmaWNhdGlvbiBmb3IgdGhpcyBkYXRhLlxu4oS5IFNwZWNpZnkgdGhlIGNvbHVtbiB0eXBlcyBvciBzZXQgYHNob3dfY29sX3R5cGVzID0gRkFMU0VgIHRvIHF1aWV0IHRoaXMgbWVzc2FnZS5cbiJ9 -->

```
Rows: 838 Columns: 36
── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Delimiter: "\t"
chr (26): Kids_First_Participant_ID, Kids_First_Biospecimen_ID_normal, sampl...
dbl (10): broad_histology_order, age_at_diagnosis_days, age_last_update_days...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



<!-- rnb-message-end -->

<!-- rnb-source-begin eyJkYXRhIjpbImNwZ3MgPC0gcmVhZF9saW5lcyhjcGdfZmlsZSkiLCIiLCJtYWYgPC0gZGF0YS50YWJsZTo6ZnJlYWQoZmlsZS5wYXRoKHJvb3RfZGlyLCBwYXJhbXMkbWFmX2luKSwgZGF0YS50YWJsZSA9IEYpICU+JVxuICBmaWx0ZXIoVHVtb3JfU2FtcGxlX0JhcmNvZGUgJWluJSB0dW1vcl9pZHMsXG4gICAgICAgICBIdWdvX1N5bWJvbCAlaW4lIGNwZ3MpICU+JVxuICB1bmlxdWUoKSAlPiVcbiAgd3JpdGVfdHN2KGZpbGUucGF0aChhbmFseXNpc19kaXIsIHBhcmFtcyRtYWZfZ29pKSkiXX0= -->

```r
cpgs <- read_lines(cpg_file)

maf <- data.table::fread(file.path(root_dir, params$maf_in), data.table = F) %>%
  filter(Tumor_Sample_Barcode %in% tumor_ids,
         Hugo_Symbol %in% cpgs) %>%
  unique() %>%
  write_tsv(file.path(analysis_dir, params$maf_goi))
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->

