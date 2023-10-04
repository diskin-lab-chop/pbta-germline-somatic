#' Create data summary
#'
#' @param x
#'
#' @return data summary table
#' @export
#'
#' @examples

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#' Create plot_exposure_violin
#'
#' @param df
#' @param x 
#' @param y
#' @param sig
#'
#' @return violin plot
#' @export
#'
#' @examples
plot_exposure_violin <- function(df, x, y, sig){
  
  comparisons_mmr <- list(c("MMR", "BRCA/interacting"), c("MMR", "Other repair"),
                          c("MMR", "Non-repair"))
  
  comparisons_brca <- list(c("BRCA/interacting", "MMR"), c("BRCA/interacting", "Other repair"), 
                           c("BRCA/interacting", "Non-repair"))
  
  vplot <- df %>%
    # generate violin plot and save 
    ggplot(aes(x = !!rlang::ensym(x), y = !!rlang::ensym(y), fill = !!rlang::ensym(x))) +
    geom_violin(binaxis = "y", stackdir = "center") +
    geom_boxplot(width=0.1) + 
    labs(x = 'Germline Variant Class', y = paste0(sig, " Exposure")) +
    theme_minimal() +
    theme(legend.position = 'none',
          text = element_text(size = 12)) +
    scale_x_discrete(labels = c(paste0("MMR\n (n=", length(mmr_ids), ")"), 
                                paste0("BRCA/\nBRCA-interacting\n (n=", length(brca_ids), ")"),
                                paste0("Other repair\n (n=", length(otherRepair_ids), ")"),
                                paste0("No DNA repair\n (n=", length(ctrl_ids), ")"))) +
    scale_fill_npg() +
    stat_summary(fun.data=data_summary,
                 show.legend = F) 
    
    if (sig == "SBS3"){
      
      vplot <- vplot + stat_compare_means(comparisons = comparisons_brca, size = 3)
      
    } else {
      
      vplot <- vplot + stat_compare_means(comparisons = comparisons_mmr, size = 3)
      
    }

  
  return(vplot)
  
}


