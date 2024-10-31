
# calculate P-LP variant enrichment in case cohort (here, PBTA) relative to control cohort
calculate_enrichment <- function(df) {
  
  # Loop through rows of df to calcualte enrichment 
  for (j in 1:nrow(df)){
    
    # run Fisher's exact test on contingency table
    fisher_test <- fisher.test(matrix(c(df$count_with_plp_case[j],
                                        df$count_without_plp_case[j],
                                        df$count_with_plp_control[j],
                                        df$count_without_plp_control[j]),
                                      2, 2))
    
    # Add Fisher's exact test results to df
    df$OR[j] = fisher_test$estimate
    df$p[j] = fisher_test$p.value
    df$ci.int1[j] = fisher_test$conf.int[1]
    df$ci.int2[j] = fisher_test$conf.int[2]
    
  }
  
  # Correct for multiple tests
  df <- df %>%
    dplyr::mutate(padj = p.adjust(p, method = "bonferroni"))
  
  # return df with enrichment calculations
  return(df)
  
}

## Plotting functions

# plot p-value function 
plot_pvalue <- function(enr_df, facet_var,
                        to_retain = NULL){
  
  # get ntests to plot bonferroni-adjusted sig threshold
  ntests <- length(unique(enr_df[[facet_var]]))
  
  # filter enrichment results for levels of variable of interest, when `to_retain` provided
  if (!is.null(to_retain)){
    
    enr_df <- enr_df %>%
      dplyr::filter(!!sym(facet_var) %in% to_retain)
    
  }
  
  pval_plot <- enr_df %>% 
    # log-transform
    mutate(p = -log10(p)) %>%
    
    ggplot(aes(x = p, y = factor(cohort))) +
    geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
    labs(x = "-log10(p)\n", y = "") + 
    # add bonferroni-adjusted signicance threshold
    geom_vline(xintercept = -log10(0.05/ntests), linetype = "dashed") + 
    xlim(0, NA) + 
    # facet wrap by variable, placing strips to the left of plots
    facet_wrap(facet_var, scale = "fixed",
               nrow = length(unique(enr_df[[facet_var]])),
               labeller = labeller(.cols = label_wrap_gen(14)),
               strip.position = "left") +
    theme_Publication() +
    # rotate strip text 
    theme(plot.margin = unit(c(2,1,1,0), "lines"),
          strip.text.y.left = element_text(angle = 0)) +
    theme(strip.placement = "outside")
  
  return(pval_plot)
  
}

# Odds Ratio plot
plot_enr <- function(enr_df, facet_var, log_scale = FALSE){
  
  # plot log10-transformed OR, when indicated 
  if (log_scale == TRUE){
    
    enr_plot <- enr_df %>% 
      dplyr::mutate(OR = case_when(
        is.infinite(OR) ~ 1000,
        TRUE ~ OR
      )) %>%
      dplyr::mutate(ci.int2 = case_when(
        is.infinite(ci.int2) ~ 10000,
        TRUE ~ ci.int2
      )) %>%
      ggplot(aes(x = log10(OR), y = factor(cohort))) +
      geom_point(size = 3, color = "#00A087FF",
                 show.legend = FALSE) + 
      # add log10-transformed confidence intervals
      geom_errorbar(aes(xmin = log10(ci.int1), xmax = log10(ci.int2)), width = 0.2, 
                    show.legend = FALSE, color = "#00A087FF") +
      labs(x = "log10-Odds Ratio\n(95% CI)", y = NULL) + 
      # Add OR threshold line
      geom_vline(xintercept = 0, linetype = "dashed") +
      # remove y axis labels 
      scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                                "PMBB" = NULL)) +
      # facet wrap by variable, place strips to left
      facet_wrap(facet_var, scale = "fixed",
                 nrow = length(unique(enr_df[[facet_var]])),
                 #    labeller = labeller(.cols = label_wrap_gen(14)),
                 strip.position = "left") +
      expand_limits(y=0) +
      theme_Publication() +
      # remove facet wrap strips
      theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      theme(strip.placement = "outside")
    
    # plot raw OR when log10 not specified 
  } else {
    
    # cap upper CI at 50
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
      labs(x = "Odds Ratio (95% CI)\n", y = NULL) + 
      # Add OR threshold line
      geom_vline(xintercept = 1, linetype = "dashed") +
      # remove y axis labels 
      scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                                "PMBB" = NULL)) +
      # facet wrap by variable, place strips to left
      facet_wrap(facet_var, scale = "fixed",
                 nrow = length(unique(enr_df[[facet_var]])),
                 #    labeller = labeller(.cols = label_wrap_gen(14)),
                 strip.position = "left") +
      expand_limits(y=0) +
      theme_Publication() +
      # remove facet wrap strips
      theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      theme(strip.placement = "outside")
  }
  
  return(enr_plot)
  
}

# plot percent P-LP by cohort
plot_perc <- function(enr_df, facet_var){
  
  # calculate position to plot fractions to the right of plot
  x_pos <- max(enr_df$perc_plp)*1.05
  
  
  perc_plot <- enr_df %>%
    ggplot(aes(x = perc_plp, y = factor(cohort), label = fraction)) +
    geom_bar(stat = "identity", color = "black",
             show.legend = TRUE, fill = "#00A087FF") + 
    # add fractions to the right of plots
    geom_text(x = x_pos, hjust = 0, size = 4, fontface = 2) +
    labs(x = "% Cohort P-LP\n", y = NULL, fill = NULL) + 
    scale_y_discrete(labels=c("PBTA" = NULL, "gnomAD" = NULL,
                              "PMBB" = NULL)) +
    guides(fill = guide_legend(nrow = 1)) +
    # facet wrap by variable, place strip to left of plot
    facet_wrap(facet_var, scale = "fixed",
               nrow = length(unique(enr_df[[facet_var]])),
               labeller = labeller(.cols = label_wrap_gen(18)),
               strip.position = "left") +
    #  expand_limits(x=3) +
    coord_cartesian(clip = 'off') +
    theme_Publication() +
    # remove facet strips
    theme(plot.margin = unit(c(2,5,1,1), "lines"),
          legend.position = c(0.5, 1.07),
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    theme(strip.placement = "outside")
  
  return(perc_plot)
  
}