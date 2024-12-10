## Function to plot TPM z-scores against alternative splicing event PSI z-scores 


plot_tpm_vs_psi <- function(df, xlab){
  
  tpm_vs_psi_plot <- ggplot(df, aes(x = sample_PSI_zscore, y = expr_zscore)) +
    geom_point() +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", 
                colour = "red",
                fill = "pink",
                linetype="dashed") +
    labs(x = xlab,
         y = "TPM z-score") + 
    stat_cor(method = "pearson",
             label.x.npc = "left",
             label.y.npc = "top",
             size = 3) +
    facet_wrap(~variant_classification, nrow = 5) + 
    ylim(c(NA, 5)) +
    theme_Publication()
  
  return(tpm_vs_psi_plot)
  
}



plot_prot_vs_psi <- function(df, xlab){
  
  prot_vs_psi_plot <- ggplot(df, aes(x = sample_PSI_zscore, y = expr_zscore_proteo)) +
    geom_point() +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", 
                colour = "red",
                fill = "pink",
                linetype="dashed") +
    labs(x = xlab,
         y = "TPM z-score") + 
    stat_cor(method = "pearson",
             label.x.npc = "left",
             label.y.npc = "top",
             size = 3) +
    facet_wrap(~variant_classification, nrow = 5) + 
    ylim(c(NA, 5)) +
    theme_Publication()
  
  return(prot_vs_psi_plot)
  
}
