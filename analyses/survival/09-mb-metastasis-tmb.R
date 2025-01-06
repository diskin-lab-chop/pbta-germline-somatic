# R Corbett 2024
#
# Generate MB metastasis logistic regression model and plot ORs

# Load packages
library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)
library(ggpubr)

# Establish directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "survival")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# Load survival functions
source(file.path(root_dir, "figures", "theme.R"))

# wrangle data

mb_hist_file <- file.path(results_dir, 
                          "mb-histologies-plus-shh-methyl-subtypes-metastsis.tsv")

mb_hist <- read_tsv(mb_hist_file)

# redefine metastasis column 
mb_hist <- mb_hist %>%
  dplyr::mutate(meta = case_when(
    metastasis == "Metastasis" ~ "Yes",
    metastasis == "No metastasis" ~ "No"
  )) %>%
  dplyr::mutate(meta = factor(meta, levels = c("No", "Yes"))) %>%
  dplyr::mutate(cpg_plp = fct_relevel(cpg_plp,
                                      c("No CPG P/LP", "CPG P/LP"))) %>%
  dplyr::mutate(molecular_subtype = fct_relevel(molecular_subtype,
                                                c("MB, WNT", "MB, Group3",
                                                  "MB, Group4", "MB, SHH"))) %>%
  dplyr::mutate(age_at_diagnosis_years = age_at_diagnosis_days/365.25)

# generate logistic regression model including appropriate covariates
model <- glm(
  meta ~ extent_of_tumor_resection + molecular_subtype + cpg_plp + age_at_diagnosis_years,
  data = mb_hist[!mb_hist$molecular_subtype %in% c("MB, To be classified"),],
  family = binomial
)
summary(model)


### Generate Forest plot

# Set up ordering and labels for y-axis
term_order <- rev(paste0(unlist(lapply(names(model$xlevels), function(x) rep(x, length(model$xlevels[[x]])))),
                         as.vector(unlist(model$xlevels))))

term_labels <- rev(as.vector(unlist(model$xlevels)))

numeric_terms <- names(model$coefficients)[!names(model$coefficients) %in% term_order]

term_order <- c(numeric_terms, term_order)

term_labels <- c(numeric_terms, term_labels)

# Convert metastasis model result to data frame, and exponentiate estimates/CIs to get ORs
survival_df <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  bind_cols(matrix(exp(confint(model)), ncol = 2, dimnames = list(NULL, c("lower .95", "upper .95")))) %>%
  dplyr::filter(!is.nan(`Pr(>|z|)`)) %>%
  dplyr::mutate(`exp(coef)` = exp(Estimate)) %>%
  
  # Add references
  add_row(term = term_order[!term_order %in% broom::tidy(model)$term], 
          `exp(coef)` = 1) %>%
  mutate(
    conf.low = `lower .95`,
    conf.high = `upper .95`,
    estimate = `exp(coef)`,
    p.value = `Pr(>|z|)`,
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  ) %>%
  filter(estimate > 1e-4 & estimate < 1500,
         term != "(Intercept)")

if (length(names(model$xlevels)) > 0) {
  
  survival_df <- survival_df %>%
    arrange(term) %>%
    dplyr::mutate(term = str_replace_all(term, paste(names(model$xlevels), collapse = "|"), "")) %>%
    dplyr::mutate(term = str_replace_all(term, "_", " ")) %>%
    dplyr::mutate(term = fct_relevel(term, unique(term)))
  
}

forest_plot <- ggplot(survival_df) +
  aes(x = estimate, y = term, fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(xmin = conf.low,xmax = conf.high,
    ), height = 0.15, linewidth = 0.65) + 
  geom_point(size = 3.5, shape = 23) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(xintercept = 1, linetype = 3
  ) +
  labs(x = "Odds Ratio ± 95% CI", y = "",
       subtitle = glue::glue("Metastasis: N = 96 with 43 events")
  ) + 
  # log-scale the x-axis
  #  scale_x_log10() +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold"),
    plot.margin = margin(r=6, unit = "pt")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()

# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.01, 
      paste0("P = ", format(round(p.value, 2), nsmall = 2)),
      "P < 0.01"
    ),
    conf.low = format(round(conf.low, 2), nsmall = 2),
    conf.high = format(round(conf.high, 2), nsmall = 2),
    estimate = format(round(estimate, 2), nsmall = 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  dplyr::mutate(hr_ci = str_replace_all(hr_ci, " - ", "-"),
                hr_ci = str_replace_all(hr_ci, "  ", ""),
                hr_ci = str_replace_all(hr_ci, "- ", "-")) %>%
  dplyr::select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = value)

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0, size = 3,
            nudge_x = -0.5) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get covered
    plot.margin = margin(6, 0, 6, 0, unit = "pt"),
    #  plot.subtitle = element_text(face = "bold")
  ) +
  scale_x_discrete(labels = c("    OR (95% CI)          ", 
                              "P-value              "),
                   position = "top")

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), 
                                    scale = 1, align = "h", hjust = 0, ncol = 2)

pdf(NULL)

pdf(file.path(plot_dir,
              "forest-MB-metastasis-additive-resection-subtype-plp-age-dx.pdf"),
    width = 8, height = 3)

forest_panels

dev.off()

# Plot TMB by subtype and CPG P/LP status

tmb <- read_tsv(file.path(data_dir, "snv-mutation-tmb-coding.tsv")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_tumor = Tumor_Sample_Barcode)

mb_hist <- mb_hist %>%
  left_join(tmb %>% dplyr::select(Kids_First_Biospecimen_ID_tumor,
                                  tmb))

mb_hist %>%
  dplyr::filter(molecular_subtype != "MB, To be classified") %>%
  
  ggplot(aes(x = cpg_plp, y = log10(tmb), fill = cpg_plp)) + 
  geom_jitter(shape = 21, 
              width = 0.2, size = 2.5, alpha = 0.85,
              show.legend = FALSE) +
  geom_boxplot(alpha = 0.05, outlier.shape = NA,
               show.legend = FALSE) +
  stat_compare_means(method = "wilcox",
                     comparisons = list(c("CPG P/LP", "No CPG P/LP")),
                     method.args = list(alternative = "less")) + 
  labs(x = NULL, y = bquote(bold(log[10](TMB~"[Mut/Mb]")))) +
  facet_wrap(~molecular_subtype, nrow = 2, scale = "free_y") +
  scale_y_continuous(expand = expansion(mult = .2)) +
  theme_Publication()

ggsave(file.path(plot_dir, "mb-tmb-by-subtype-cpg-plp-status.pdf"),
       width = 6, height = 5)

# Plot age at dx by subtype and CPG P/LP status

mb_hist %>%
  dplyr::filter(molecular_subtype != "MB, To be classified") %>%
  
  ggplot(aes(x = cpg_plp, y = age_at_diagnosis_years, fill = cpg_plp)) + 
  geom_jitter(shape = 21, 
              width = 0.2, size = 2.5, alpha = 0.85,
              show.legend = FALSE) +
  geom_boxplot(alpha = 0.05, outlier.shape = NA,
               show.legend = FALSE) +
  stat_compare_means(method = "wilcox",
                     comparisons = list(c("CPG P/LP", "No CPG P/LP")),
                     method.args = list(alternative = "two.sided")) + 
  labs(x = NULL, y = "Age at diagnosis (years)") +
  facet_wrap(~molecular_subtype, nrow = 2, scale = "free_y") +
  scale_y_continuous(expand = expansion(mult = .2)) +
  theme_Publication()

ggsave(file.path(plot_dir, "mb-age-dx-by-subtype-cpg-plp-status.pdf"),
       width = 6, height = 5)

