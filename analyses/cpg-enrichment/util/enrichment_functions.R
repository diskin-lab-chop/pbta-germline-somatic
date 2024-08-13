


plot_pvalue <- function(enr_df, facet_var){
  
  ntests <- length(unique(enr_df[[facet_var]]))
  
  pval_plot <- enr_df %>% 
    mutate(p = -log10(p)) %>%
    ggplot(aes(x = p, y = factor(cohort))) +
    geom_point(size = 3, show.legend = FALSE, color = "#00A087FF") + 
    labs(x = "-log10(p)", y = "") + 
    geom_vline(xintercept = -log10(0.05/ntests), linetype = "dashed") + 
    facet_wrap(facet_var, scale = "fixed",
               nrow = length(unique(enr_df[[facet_var]])),
               labeller = labeller(.cols = label_wrap_gen(18))) +  
    theme_Publication() +
    theme(plot.margin = unit(c(2,1,1,0), "lines")) +
    theme(strip.placement = "outside")
  
  return(pval_plot)
  
}


plot_enr <- function(enr_df, facet_var, log_scale = FALSE){
  
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
               labeller = labeller(.cols = label_wrap_gen(18))) +  
    expand_limits(y=0) +
    theme_Publication() +
    theme(plot.margin = unit(c(2,0.5,1,0.25), "lines")) +
    theme(strip.placement = "outside")
  
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
               labeller = labeller(.cols = label_wrap_gen(18))) +
    expand_limits(x=3) +
    coord_cartesian(clip = 'off') +
    theme_Publication() +
    theme(plot.margin = unit(c(2,5,1,1), "lines"),
          legend.position = c(0.5, 1.07)) +
    theme(strip.placement = "outside")
  
  return(perc_plot)
  
}