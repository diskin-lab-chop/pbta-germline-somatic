# Add CNVkit results to df by sample ID and gene ID

add_cnvkit <- function(loh_df, cnvkit_df){
  
  # add empty columns to loh_df
  loh_df <- loh_df %>%
    dplyr::mutate(cnvkit_log2 = NA_integer_,
                  cnvkit_cn = NA_integer_,
                  cnvkit_cn1 = NA_integer_,
                  cnvkit_cn2 = NA_integer_)
  
  # loop through loh_df rows to add cnvkit data
  for (i in 1:nrow(loh_df)){
    
    # extract sample and variant coordinates
    id <- loh_df$Kids_First_Biospecimen_ID_tumor[i]
    variant_chr <- loh_df$chr[i]
    variant_start <- loh_df$start[i]

    # subset cnvkit for sample and position
    cnvkit_filtered <- cnvkit_df %>%
      dplyr::filter(ID == id,
                    chromosome == variant_chr,
                    variant_start > start & variant_start < end)
    
    # add cnvkit data
    if (nrow(cnvkit_filtered) > 0){
      
      loh_df$cnvkit_log2[i] <- cnvkit_filtered$log2
      loh_df$cnvkit_cn[i] <- cnvkit_filtered$cn
      loh_df$cnvkit_cn1[i] <- cnvkit_filtered$cn1
      loh_df$cnvkit_cn2[i] <- cnvkit_filtered$cn2
      
    }
    
  }
  
  return(loh_df)
  
}