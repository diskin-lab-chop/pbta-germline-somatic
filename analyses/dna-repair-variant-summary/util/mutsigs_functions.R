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
plot_exposure_violin <- function(df, x, y, sig, group){
  
  comparisons <- ifelse(sig == "SBS3",
                        list(c("HR", "Non-HR")),
                        ifelse(sig == "SBS30",
                               list(c("BER", "Non-BER")),
                               list(c("MMR", "Non-MMR"))))
  
  brca_labels <- c(paste0("HR P/LP\n (n=", sum(df$germline_variant == "HR"), ")"), 
                   paste0("HR WT\n (n=", sum(df$germline_variant == "Non-HR"), ")"))
  
  mmr_labels <- c(paste0("MMR P/LP\n (n=", sum(df$germline_variant == "MMR"), ")"), 
                  paste0("MMR WT\n (n=", sum(df$germline_variant == "Non-MMR"), ")"))
  
  ber_labels <- c(paste0("BER P/LP\n (n=", sum(df$germline_variant == "BER"), ")"), 
                  paste0("BER WT\n (n=", sum(df$germline_variant == "Non-BER"), ")"))
  
  x_labels <- if (sig == "SBS3"){
    
    brca_labels
    
  } else if (sig == "SBS30"){
    
    ber_labels
    
  } else {
    
    mmr_labels
    
  }
  
  vplot <- df %>%
    # generate violin plot and save 
    ggplot(aes(x = !!rlang::ensym(x), y = !!rlang::ensym(y), fill = !!rlang::ensym(x))) +
    geom_violin(binaxis = "y", stackdir = "center", 
                show.legend = FALSE) +
    geom_boxplot(width=0.1, show.legend = FALSE) + 
    labs(x = 'Germline Variant Class', y = paste0(sig, " Exposure"),
         title = group) +
    #  theme_minimal() +
    # theme(legend.position = 'none',
    #       text = element_text(size = 12),
    #       plot.title = element_text(size = 8, hjust = 0.5)) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_npg() +
    stat_summary(fun.data=data_summary,
                 show.legend = F) +
    scale_y_continuous(expand = expansion(mult = .2)) +
    theme_Publication() +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          axis.title.x = element_text(size = 13)) 
  
  vplot <- vplot + stat_compare_means(method = "wilcox",
                                      comparisons = comparisons, size = 3,
                                      method.args = list(alternative = "greater"))
  
  return(vplot)
  
}


