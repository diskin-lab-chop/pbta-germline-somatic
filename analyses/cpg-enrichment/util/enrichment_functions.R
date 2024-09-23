
calculate_enrichment <- function(df) {
  
  for (j in 1:nrow(df)){
    
    fisher_test <- fisher.test(matrix(c(df$count_with_plp_case[j],
                                        df$count_without_plp_case[j],
                                        df$count_with_plp_control[j],
                                        df$count_without_plp_control[j]),
                                      2, 2))
    
    df$OR[j] = fisher_test$estimate
    df$p[j] = fisher_test$p.value
    df$ci.int1[j] = fisher_test$conf.int[1]
    df$ci.int2[j] = fisher_test$conf.int[2]
    
  }
  
  # Correct for multiple tests
  df <- df %>%
    dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))
  
  return(df)
  
}


plot_pvalue <- function(enr_df, facet_var,
                        to_retain = NULL){
  
  ntests <- length(unique(enr_df[[facet_var]]))
  
  if (!is.null(to_retain)){
    
    enr_df <- enr_df %>%
      dplyr::filter(!!sym(facet_var) %in% to_retain)
    
  }
  
  pval_plot <- enr_df %>% 
    mutate(p = -log10(p)) %>%
    ggplot(aes(x = p, y = factor(cohort))) +
    geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
    labs(x = "-log10(p)", y = "") + 
    geom_vline(xintercept = -log10(0.05/ntests), linetype = "dashed") + 
    xlim(0, NA) + 
    facet_wrap(facet_var, scale = "fixed",
               nrow = length(unique(enr_df[[facet_var]])),
               labeller = labeller(.cols = label_wrap_gen(14)),
               strip.position = "left") +
    theme_Publication() +
    theme(plot.margin = unit(c(2,1,1,0), "lines"),
          strip.text.y.left = element_text(angle = 0)) +
    theme(strip.placement = "outside")
  
  return(pval_plot)
  
}


plot_enr <- function(enr_df, facet_var, log_scale = FALSE){

  if (log_scale == TRUE){
    
    # enr_df <- enr_df %>%
    #   dplyr::mutate(ci.int2 = case_when(
    #     ci.int2 > 1000 ~ 1000,
    #     TRUE ~ ci.int2
    #   ))
  
    enr_plot <- enr_df %>% 
      ggplot(aes(x = log10(OR), y = factor(cohort))) +
      geom_point(size = 3, color = "#00A087FF",
                 show.legend = FALSE) + 
      geom_errorbar(aes(xmin = log10(ci.int1), xmax = log10(ci.int2)), width = 0.2, 
                    show.legend = FALSE, color = "#00A087FF") +
      labs(x = "log10-Odds Ratio (95% CI)", y = NULL) + 
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_y_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                                "PMBB" = "")) +
      facet_wrap(facet_var, scale = "fixed",
                 nrow = length(unique(enr_df[[facet_var]])),
                 labeller = labeller(.cols = label_wrap_gen(18)),
                 strip.position = "left") +
      expand_limits(y=0) +
      theme_Publication() +
      theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      theme(strip.placement = "outside")
  
  } else {
    
    enr_df <- enr_df %>%
      dplyr::mutate(ci.int2 = case_when(
        ci.int2 > 50 ~ 50,
        TRUE ~ ci.int2
      ))
    
    enr_plot <- enr_df %>% 
      ggplot(aes(x = OR, y = factor(cohort))) +
      geom_point(size = 3, color = "#00A087FF",
                 show.legend = FALSE) + 
      geom_errorbar(aes(xmin = ci.int1, xmax = ci.int2), width = 0.2, 
                    show.legend = FALSE, color = "#00A087FF") +
      labs(x = "Odds Ratio (95% CI)", y = NULL) + 
      geom_vline(xintercept = 1, linetype = "dashed") +
      scale_y_discrete(labels=c("PBTA" = "", "gnomAD" = "",
                                "PMBB" = "")) +
      facet_wrap(facet_var, scale = "fixed",
                 nrow = length(unique(enr_df[[facet_var]])),
                 labeller = labeller(.cols = label_wrap_gen(18)),
                 strip.position = "left") +
      expand_limits(y=0) +
      theme_Publication() +
      theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      theme(strip.placement = "outside")
  }
  
  return(enr_plot)
  
}


plot_perc <- function(enr_df, facet_var){
  
  x_pos <- max(enr_df$perc_plp)*1.05
  
  perc_plot <- enr_df %>%
    ggplot(aes(x = perc_plp, y = factor(cohort), label = fraction)) +
    geom_bar(stat = "identity", color = "black",
             show.legend = TRUE, fill = "#00A087FF") + 
    geom_text(x = x_pos, hjust = 0, size = 4, fontface = 2) +
    labs(x = "% Cohort P/LP", y = NULL, fill = NULL) + 
    scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                              "PMBB" = NULL)) +
    guides(fill = guide_legend(nrow = 1)) +
    facet_wrap(facet_var, scale = "fixed",
               nrow = length(unique(enr_df[[facet_var]])),
               labeller = labeller(.cols = label_wrap_gen(18)),
               strip.position = "left") +
  #  expand_limits(x=3) +
    coord_cartesian(clip = 'off') +
    theme_Publication() +
    theme(plot.margin = unit(c(2,5,1,1), "lines"),
          legend.position = c(0.5, 1.07),
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    theme(strip.placement = "outside")
  
  return(perc_plot)
  
}