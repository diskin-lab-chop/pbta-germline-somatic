## GSVA plotting functions

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

#' Create gsva score difference barplot between two groups
#'
#' @param gsva
#' @param x 
#' @param y
#' @param plot_name
#'
#' @return barplot plot
#' @export
#'
#' @examples

gsva_barplot <- function(gsva, x, y, plot_name){
  
  bplot <- gsva %>%
    ggplot(aes(x = !!rlang::ensym(x),
               y = reorder(!!rlang::ensym(y), -!!rlang::ensym(x)),
               label = ifelse(pvalue < 0.05, "*", ""))) +
    geom_bar(stat="identity", color = "black") +
    labs(x = 'GSEA score difference', y = NULL,
         title = glue::glue("HGG, {plot_name} \n germline P/LP vs. no P/LP")) +
    geom_vline(xintercept = 0) + 
    xlim(-0.4, 0.4) +
    geom_text(hjust = -0.5, vjust = 0.75, size = 7) +
    scale_fill_npg() + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(bplot)
  
}

#' Create gsva score violin plot
#'
#' @param gsva
#' @param x 
#' @param y
#' @param title
#'
#' @return violin plot
#' @export
#'
#' @examples

gsva_violin <- function(gsva, x, y, title){
  
  comparisons <- list(c("MMR", "BRCA/BRCA-interacting"), c("MMR", "Other repair"),
                      c("MMR", "Non-repair"))
  
  vplot <- gsva %>%
    ggplot(aes(x = !!rlang::ensym(x), y = !!rlang::ensym(y), fill = !!rlang::ensym(x))) +
    geom_violin(binaxis = "y", stackdir = "center",
                show.legend = FALSE) +
    geom_boxplot(width=0.1,
                 show.legend = FALSE) + 
    labs(x = 'Germline Variant Class', y = " GSVA score",
         title = title) +
    scale_x_discrete(labels = c(paste0("MMR\n (n=", length(unique(mmr_ids)), ")"), 
                                paste0("BRCA/\nBRCA-interacting\n (n=", length(unique(brca_ids)), ")"), 
                                paste0("Other repair\n (n=", length(unique(otherRepair_ids)), ")"),
                                paste0("No DNA repair\n (n=", length(unique(ctrl_ids)), ")"))) +
    scale_fill_npg() +
    stat_summary(fun.data=data_summary,
                 show.legend = F) + 
    stat_compare_means(comparisons = comparisons, size = 3) +
    theme_Publication()
  
  
}
