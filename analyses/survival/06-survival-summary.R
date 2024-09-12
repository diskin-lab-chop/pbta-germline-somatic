# R Corbett 2023
#
# Generate survival summary table for PBTA Ancestry project

# Load packages
library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)

# Establish directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival")
input_dir <- file.path(analysis_dir, "results")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# Load survival functions
source(file.path(analysis_dir, "util", "survival_models.R"))
source(file.path(root_dir, "figures", "theme.R"))

hist_file <- file.path(input_dir, 
                           "germline-primary-plus-tumor-histologies-plot-groups-clin-meta-subtype.tsv")

# Define df of histology names and folder names
hist_df <- data.frame(group = c("Atypical Teratoid Rhabdoid Tumor", "DIPG or DMG", 
                                "High-grade glioma", 
                                "Ependymoma", "Mixed neuronal-glial tumor", 
                                "Low-grade glioma", "Medulloblastoma", "Meningioma", 
                                "Neurofibroma plexiform"),
                      hist = c("ATRT", "DMG", "HGG", "EPN", 
                               "GNG-GNT", "LGG", "MB", "MNG", "NFP"))

# Load ancestry file
hist <- read_tsv(hist_file)

# Create subtype df of molecular subgroup names and folder names
subtype_df <- hist %>%
  count(mol_sub_group, plot_group) %>%
  filter(n >=15 & !grepl("classified", mol_sub_group) & !is.na(mol_sub_group)) %>%
  dplyr::mutate(hist = unlist(lapply(strsplit(mol_sub_group, ", "), function(x) x[1])),
                subtype = case_when(
                  grepl(",", mol_sub_group) ~ unlist(lapply(strsplit(mol_sub_group, ", "), function(x) x[2])),
                  TRUE ~ NA_character_)) %>%
  dplyr::mutate(subtype = str_replace_all(subtype, " ", "-")) %>%
  dplyr::mutate(hist = case_when(
    hist == "GNG/GNT" ~ "GNG-GNT",
    TRUE ~ hist
  ))

# merge histology and subtype names 
group_df <- subtype_df %>%
  dplyr::rename(group = mol_sub_group) %>%
  dplyr::mutate(subtype = glue::glue("{hist}_{subtype}")) %>%
  dplyr::select(group, hist, subtype) %>%
  bind_rows(hist_df) %>%
  dplyr::mutate(subtype = case_when(
    is.na(subtype) ~ hist,
    TRUE ~ subtype
  )) %>%
  arrange(group)

# Create empty df to store survival model summary stats
os_stats <- data.frame(group = c(group_df$group),
                       type = rep("Overall survival", length(group_df$group)),
                       HR = rep(0, length(group_df$group)),
                       p = rep(0, length(group_df$group)),
                       CI_lower = rep(0, length(group_df$group)),
                       CI_upper = rep(0, length(group_df$group)))
                       
efs_stats <- data.frame(group = c(group_df$group),
                        type = rep("Event-free survival", length(group_df$group)),
                        HR = rep(0, length(group_df$group)),
                        p = rep(0, length(group_df$group)),
                        CI_lower = rep(0, length(group_df$group)),
                        CI_upper = rep(0, length(group_df$group)))

# Loop through histologies and molecular subtypes and extract OS and EFS hazard ratios and associated p-values for CPG P-LP carriers vs. non-carriers
for (i in 1:nrow(group_df)){
  
  input_dir <- file.path(analysis_dir, "results", group_df$hist[i])
  
  if (grepl("Low-grade glioma|LGG", group_df$group[i])){
    survival_os_plp <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    survival_os_plp <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  # tidy survival results 
  os_df <- broom::tidy(survival_os_plp)
  os_ci_df <- summary(survival_os_plp)$conf.int
  
  # extract OS HRs and p-values
  os_stats[os_stats$group == group_df$group[i],]$HR <- exp(os_df$estimate[os_df$term == "cpgPLP_statuscpgPLP"])
  os_stats[os_stats$group == group_df$group[i],]$p <- os_df$p.value[os_df$term == "cpgPLP_statuscpgPLP"]
  
  os_stats[os_stats$group == group_df$group[i],]$CI_lower <- os_ci_df["cpgPLP_statuscpgPLP", "lower .95"]
  os_stats[os_stats$group == group_df$group[i],]$CI_upper <- os_ci_df["cpgPLP_statuscpgPLP", "upper .95"]
  
  
  if (grepl("Low-grade glioma|LGG", group_df$group[i])){
    survival_efs_plp <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    survival_efs_plp <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  efs_df <- broom::tidy(survival_efs_plp)
  efs_ci_df <- summary(survival_efs_plp)$conf.int
  
  # extract EFS HRs and p-values
  efs_stats[efs_stats$group == group_df$group[i],]$HR <- exp(efs_df$estimate[efs_df$term == "cpgPLP_statuscpgPLP"])
  efs_stats[efs_stats$group == group_df$group[i],]$p <- efs_df$p.value[efs_df$term == "cpgPLP_statuscpgPLP"]
  
  efs_stats[efs_stats$group == group_df$group[i],]$CI_lower <- efs_ci_df["cpgPLP_statuscpgPLP", "lower .95"]
  efs_stats[efs_stats$group == group_df$group[i],]$CI_upper <- efs_ci_df["cpgPLP_statuscpgPLP", "upper .95"]
  
}

# Some models have massively inflated p-values driven by single events; set these to NA
os_stats[abs(log10(os_stats$HR)) > 3 & !is.na(os_stats$HR), c("HR", "p", "CI_lower", "CI_upper")] <- NA
efs_stats[abs(log10(efs_stats$HR)) > 3 & !is.na(efs_stats$HR), c("HR", "p", "CI_lower", "CI_upper")] <- NA

# merge OS and EFS stats
survival_stats <- os_stats %>%
  bind_rows(efs_stats) %>%
  dplyr::filter(!is.na(HR))

# add p-value lable column for HRs with p<0.1
survival_stats <- survival_stats %>%
  dplyr::mutate(p_label = case_when(
    p< 0.01 ~ "*p<0.01",
    p < 0.05 ~ glue::glue("*p={round(p, 2)}"),
    p < 0.1 ~ glue::glue("^p={round(p, 2)}"),
    TRUE ~ ""
  )) %>%
  dplyr::mutate(p_label = case_when(
    p< 0.01 ~ "*p<0.01",
    p < 0.05 ~ glue::glue("*p={round(p, 2)}"),
    p < 0.1 ~ glue::glue("^p={round(p, 2)}"),
    TRUE ~ ""
  ))

# reorder rows
group_order <- c("Atypical Teratoid Rhabdoid Tumor",
                 "DIPG or DMG",
                 "Ependymoma",
                 "EPN, PF A",
                 "EPN, ST ZFTA",
                 "High-grade glioma",
                 "HGG, H3 WT",
                 "Low-grade glioma",
                 "LGG, BRAF V600E",
                 "LGG, BRAF fusion",
                 "LGG, Other alteration",
                 "Medulloblastoma",
                 "MB, Group3", 
                 "MB, Group4",
                 "MB, SHH",
                 "Meningioma",
                 "Mixed neuronal-glial tumor",
                 "GNG/GNT, BRAF V600E",
                 "GNG/GNT, Other alteration",
                 "Neurofibroma plexiform")

survival_stats <- survival_stats %>% 
  dplyr::filter(group %in% group_order) %>%
  dplyr::mutate(group = fct_relevel(group, 
                                    rev(group_order)))


# Plot hazard ratios
surv_plot <- survival_stats %>%
  dplyr::mutate(CI_upper = case_when(
    CI_upper > 100 ~ 100,
    TRUE ~ CI_upper
  )) %>%
  
  ggplot(aes(x = log10(HR), y = group,
             label = p_label)) +
  geom_point(size = 3, color = "#00A087FF",
             show.legend = FALSE) + 
  geom_errorbar(aes(xmin = log10(CI_lower), xmax = log10(CI_upper)), width = 0.2, 
                show.legend = FALSE, color = "#00A087FF") +
  geom_text(x = 1, hjust = 0, size = 3.5, fontface = 2) +
  labs(x = "log10-P/LP carrier hazard ratio (95% CI)", y = NULL) + 
  xlim(-2, 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~type, nrow = 1) +
  theme_Publication() +
  theme(axis.text.x = ggtext::element_markdown())

surv_plot

# save plot
ggsave(file.path(plot_dir, "survival-hr-plp-vs-wt.pdf"),
       width = 8, height = 8)

# Save table to output
survival_stats %>% 
  dplyr::select(-p_label) %>%
write_tsv(file.path(results_dir, 
                    "median-survival-by-ancestry-cancer-group.tsv"))  

