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
    
    # extract sample and gene IDs
    id <- loh_df$Kids_First_Biospecimen_ID_tumor[i]
    gene_symbol <- loh_df$Hugo_Symbol[i]
    
    # subset cnvkit for sample and gene
    cnvkit_filtered <- cnvkit_df %>%
      dplyr::filter(ID == id,
                    grepl(glue::glue(",{gene_symbol},"), gene)) %>%
      distinct()
    
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