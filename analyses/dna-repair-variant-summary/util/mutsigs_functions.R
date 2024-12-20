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
  
  comparisons <- ifelse(sig == "SBS3",
                            list(c("HR", "Non-HR")),
                            list(c("MMR", "Non-MMR")))
  
  mmr_labels <- c(paste0("HR P/LP\n (n=", sum(df$germline_variant == "HR"), ")"), 
                  paste0("HR WT\n (n=", sum(df$germline_variant == "Non-HR"), ")"))
  
  brca_labels <- c(paste0("MMR P/LP\n (n=", sum(df$germline_variant == "MMR"), ")"), 
                   paste0("MMR WT\n (n=", sum(df$germline_variant == "Non-MMR"), ")"))
  
  x_labels <- if (sig == "SBS3"){
    
                   mmr_labels
    
  } else {
    
    brca_labels
    
  }
                   brca_labels
  
  vplot <- df %>%
    # generate violin plot and save 
    ggplot(aes(x = !!rlang::ensym(x), y = !!rlang::ensym(y), fill = !!rlang::ensym(x))) +
    geom_violin(binaxis = "y", stackdir = "center", 
                show.legend = FALSE) +
    geom_boxplot(width=0.1, show.legend = FALSE) + 
    labs(x = 'Germline Variant Class', y = paste0(sig, " Exposure"),
         title = "HGG, H3-WT") +
    theme_minimal() +
    theme(legend.position = 'none',
          text = element_text(size = 12)) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_npg() +
    stat_summary(fun.data=data_summary,
                 show.legend = F) +
    scale_y_continuous(expand = expansion(mult = .2)) +
    theme_Publication()

    vplot <- vplot + stat_compare_means(method = "wilcox",
                                        comparisons = comparisons, size = 3,
                                        method.args = list(alternative = "greater"))

  return(vplot)
  
}


