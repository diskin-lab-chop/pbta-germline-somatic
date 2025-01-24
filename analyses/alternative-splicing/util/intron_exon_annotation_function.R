library(biomaRt)

get_exon_intron_numbers <- function(df){
  
  # load ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # set start coordinates to 1-based
  df <- plp_splice_candidates %>%
    dplyr::mutate(region_start = region_start + 1,
                  upstream_exon_start = upstream_exon_start + 1,
                  downstream_exon_start = downstream_exon_start + 1)
  
  # extract gene symbols
  genes <- df %>%
    distinct(Hugo_Symbol) %>%
    pull(Hugo_Symbol)
  
  # GBA needs to be updated to GBA1
  genes <- str_replace(genes, "GBA", "GBA1")
  
  # Generate data frame of exon ids for all transcripts
  exons_feature_bm <- getBM(
    attributes = c("external_gene_name", "ensembl_transcript_id", "transcript_is_canonical", "ensembl_exon_id", "chromosome_name"),
    filters = "external_gene_name",
    values = genes, 
    mart = ensembl
  )
  
  # create separate vectors for canonical and non-canonical transcripts
  canonical_ids <- exons_feature_bm %>%
    dplyr::filter(transcript_is_canonical == 1) %>%
    distinct(ensembl_transcript_id) %>%
    pull(ensembl_transcript_id)
  
  noncanonical_ids <- exons_feature_bm %>%
    dplyr::filter(is.na(transcript_is_canonical)) %>%
    distinct(ensembl_transcript_id) %>%
    pull(ensembl_transcript_id)
  
  # Generate data frame of exon coordinates for all canonical transcripts
  canonical_exons_structure_bm <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_exon_id", "rank", "chromosome_name", "exon_chrom_start", "exon_chrom_end"),
    filters = "ensembl_transcript_id",
    values = canonical_ids,  # Chromosome name (adjust as needed)
    mart = ensembl
  )
  
  # merge feature and structure info
  canonical_exons_bm <- exons_feature_bm %>%
    dplyr::filter(transcript_is_canonical == 1) %>%
    left_join(canonical_exons_structure_bm) %>%
    dplyr::mutate(chromosome_name = glue::glue("chr{chromosome_name}"))
  
  # repeat for non-canonical transcripts
  noncanonical_exons_structure_bm <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_exon_id", "rank", "chromosome_name", "exon_chrom_start", "exon_chrom_end"),
    filters = "ensembl_transcript_id",
    values = noncanonical_ids,  # Chromosome name (adjust as needed)
    mart = ensembl
  ) %>%
    arrange(ensembl_transcript_id) %>%
    distinct(exon_chrom_end, .keep_all = TRUE)
  
  noncanonical_exons_bm <- exons_feature_bm %>%
    left_join(noncanonical_exons_structure_bm) %>%
    dplyr::filter(!is.na(rank)) %>%
    dplyr::mutate(chromosome_name = glue::glue("chr{chromosome_name}"))
  
  # first, get exon ranks for canonical transcripts for all splice juntions
  df <- df %>%
    left_join(canonical_exons_bm %>%
                dplyr::select(chromosome_name,
                              ensembl_transcript_id,
                              rank, exon_chrom_start),
              by = c("chr" = "chromosome_name",
                     "region_start" = "exon_chrom_start")) %>%
    dplyr::rename(exon_rank_canonical = rank,
                  canonical_transcript_id = ensembl_transcript_id) %>%
    left_join(canonical_exons_bm %>%
                dplyr::select(chromosome_name,
                              ensembl_transcript_id,
                              rank, exon_chrom_start),
              by = c("chr" = "chromosome_name",
                     "upstream_exon_start" = "exon_chrom_start",
                     "canonical_transcript_id" = "ensembl_transcript_id")) %>%
    dplyr::rename(upstream_exon_rank_canonical = rank) %>%
    left_join(canonical_exons_bm %>%
                dplyr::select(chromosome_name,
                              ensembl_transcript_id,
                              rank, exon_chrom_start),
              by = c("chr" = "chromosome_name",
                     "downstream_exon_start" = "exon_chrom_start",
                     "canonical_transcript_id" = "ensembl_transcript_id")) %>%
    dplyr::rename(downstream_exon_rank_canonical = rank) %>%
    # define exon and intron ranks
    dplyr::mutate(exon_intron_id = case_when(
      # if non-RI event, take exon rank
      splicing_case %in% c("SE", "A3SS", "A5SS") & !is.na(exon_rank_canonical) ~ glue::glue("exon {exon_rank_canonical}"),
      # if RI event involving single intron, intron rank will match that of upstream exon 
      splicing_case == "RI" & abs(upstream_exon_rank_canonical - downstream_exon_rank_canonical) == 1 & strand == "-" ~ glue::glue("intron {downstream_exon_rank_canonical}"),
      splicing_case == "RI" & abs(upstream_exon_rank_canonical - downstream_exon_rank_canonical) == 1 & strand == "+" ~ glue::glue("intron {upstream_exon_rank_canonical}"),
      # if RI event involves multiple introns, intron ranks will range from upstream exon rank to downstream exon rank minus one
      splicing_case == "RI" & !is.na(upstream_exon_rank_canonical) & strand == "-" ~ glue::glue("intron {downstream_exon_rank_canonical}-{upstream_exon_rank_canonical - 1}"),
      splicing_case == "RI" & !is.na(upstream_exon_rank_canonical) & strand == "+" ~ glue::glue("intron {upstream_exon_rank_canonical}-{downstream_exon_rank_canonical - 1}"),
      TRUE ~ NA
    )) %>%
      dplyr::select(-exon_rank_canonical,
                    -upstream_exon_rank_canonical,
                    -downstream_exon_rank_canonical)
  
  # for splice junctions that were not found in canonical transcripts, repeat steps above for non-canonical transcripts
  df <- df %>%
    left_join(noncanonical_exons_bm %>%
                dplyr::select(external_gene_name,
                              ensembl_transcript_id,
                              rank, exon_chrom_end),
              by = c("Hugo_Symbol" = "external_gene_name",
                     "region_end" = "exon_chrom_end")) %>%
    dplyr::rename(exon_rank_noncanonical = rank,
                  noncanonical_transcript_id = ensembl_transcript_id) %>%
    left_join(noncanonical_exons_bm %>%
                dplyr::select(external_gene_name,
                              ensembl_transcript_id,
                              rank, exon_chrom_end),
              by = c("Hugo_Symbol" = "external_gene_name",
                     "upstream_exon_end" = "exon_chrom_end",
                     "noncanonical_transcript_id" = "ensembl_transcript_id")) %>%
    dplyr::rename(upstream_exon_rank_noncanonical = rank) %>%
    left_join(noncanonical_exons_bm %>%
                dplyr::select(external_gene_name,
                              ensembl_transcript_id,
                              rank, exon_chrom_end),
              by = c("Hugo_Symbol" = "external_gene_name",
                     "downstream_exon_end" = "exon_chrom_end",
                     "noncanonical_transcript_id" = "ensembl_transcript_id")) %>%
    dplyr::rename(downstream_exon_rank_noncanonical = rank) %>%
    dplyr::mutate(exon_intron_id = case_when(
      is.na(exon_intron_id) & splicing_case %in% c("SE", "A3SS", "A5SS") & !is.na(exon_rank_noncanonical) ~ glue::glue("exon {exon_rank_noncanonical}"),
      is.na(exon_intron_id) & splicing_case == "RI" & abs(upstream_exon_rank_noncanonical - downstream_exon_rank_noncanonical) == 1 & strand == "-" ~ glue::glue("intron {downstream_exon_rank_noncanonical}"),
      is.na(exon_intron_id) & splicing_case == "RI" & abs(upstream_exon_rank_noncanonical - downstream_exon_rank_noncanonical) == 1 & strand == "+" ~ glue::glue("intron {upstream_exon_rank_noncanonical}"),
      is.na(exon_intron_id) & splicing_case == "RI" & strand == "-" ~ glue::glue("intron {downstream_exon_rank_noncanonical}-{upstream_exon_rank_noncanonical - 1}"),
      is.na(exon_intron_id) & splicing_case == "RI" & strand == "+" ~ glue::glue("intron {upstream_exon_rank_noncanonical}-{downstream_exon_rank_noncanonical - 1}"),
      TRUE ~ exon_intron_id
    )) %>%
    dplyr::mutate(ensembl_transcript_id = case_when(
      !is.na(canonical_transcript_id) ~ canonical_transcript_id,
      TRUE ~ noncanonical_transcript_id
    )) %>%
    dplyr::select(-exon_rank_noncanonical,
                  -upstream_exon_rank_noncanonical,
                  -downstream_exon_rank_noncanonical,
                  -noncanonical_transcript_id,
                  -canonical_transcript_id)
  
  return(df)
  
}
  



