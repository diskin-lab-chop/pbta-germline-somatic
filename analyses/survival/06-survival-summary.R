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
                                "Low-grade glioma", "Medulloblastoma", "Meningioma"),
                      hist = c("ATRT", "DMG", "HGG", "EPN", 
                               "GNG-GNT", "LGG", "MB", "MNG"))

# Load ancestry file
hist <- read_tsv(hist_file)

# Create subtype df of molecular subgroup names and folder names
subtype_df <- hist %>%
  # only retain patients with survival data
  dplyr::filter(!is.na(EFS_days) | !is.na(OS_days),
                mol_sub_group != "GNG/GNT, WT") %>%
  dplyr::mutate(mol_sub_group = case_when(
    grepl("DMG|HGG, H3", molecular_subtype) ~ molecular_subtype,
    TRUE ~ mol_sub_group
  )) %>%
  dplyr::count(mol_sub_group, plot_group) %>%
  # filter for groups with >=15 pts
  filter(n >=8 & !grepl("classified", mol_sub_group) & !is.na(mol_sub_group)) %>%
  dplyr::mutate(hist = unlist(lapply(strsplit(mol_sub_group, ", "), function(x) x[1])),
                subtype = str_remove(mol_sub_group, "^[^ ]+ ")
  ) %>%
  dplyr::mutate(subtype = str_replace_all(subtype, " |, ", "-")) %>%
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
                       group_n = NA_integer_,
                       type = rep("Overall survival", length(group_df$group)),
                       median_surv_years_plp = NA_integer_,
                       median_surv_years_no_plp = NA_integer_, 
                       km_p = NA_integer_,
                       HR = NA_integer_,
                       p = NA_integer_,
                       CI_lower = NA_integer_,
                       CI_upper = NA_integer_)
                       
efs_stats <- data.frame(group = c(group_df$group),
                        group_n = NA_integer_,
                        type = rep("Event-free survival", length(group_df$group)),
                        median_surv_years_plp = NA_integer_,
                        median_surv_years_no_plp = NA_integer_,
                        km_p = NA_integer_,
                        HR = NA_integer_,
                        p = NA_integer_,
                        CI_lower = NA_integer_,
                        CI_upper = NA_integer_)

# Loop through histologies and molecular subtypes and extract OS and EFS hazard ratios and associated p-values for CPG P-LP carriers vs. non-carriers
for (i in 1:nrow(group_df)){
  
  input_dir <- file.path(analysis_dir, "results", group_df$hist[i])
  
  km_os_plp <- read_rds(file.path(input_dir,
                                  glue::glue("logrank_{group_df$subtype[i]}_OS_cpgPLPstatus.RDS")))
  
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
  
  os_n <- sum(summary(km_os_plp$model)$table[,"records"])
  
  os_stats[os_stats$group == group_df$group[i],]$median_surv_years_plp <- summary(km_os_plp$model)$table[,"median"][2]/365.25
  os_stats[os_stats$group == group_df$group[i],]$median_surv_years_no_plp <- summary(km_os_plp$model)$table[,"median"][1]/365.25
  
  os_diff_obj <- survdiff(survival::Surv(OS_days, OS_status) ~ cpgPLP_status,  
                                 km_os_plp$original_data)
  
  os_stats[os_stats$group == group_df$group[i],]$km_p <- 1 - pchisq(os_diff_obj$chisq, length(os_diff_obj$n) - 1)
  
  # extract OS HRs and p-values
  os_stats[os_stats$group == group_df$group[i],]$HR <- exp(os_df$estimate[os_df$term == "cpgPLP_statusCPG P/LP"])
  os_stats[os_stats$group == group_df$group[i],]$p <- os_df$p.value[os_df$term == "cpgPLP_statusCPG P/LP"]
  
  os_stats[os_stats$group == group_df$group[i],]$CI_lower <- os_ci_df["cpgPLP_statusCPG P/LP", "lower .95"]
  os_stats[os_stats$group == group_df$group[i],]$CI_upper <- os_ci_df["cpgPLP_statusCPG P/LP", "upper .95"]
  
  
  km_efs_plp <- read_rds(file.path(input_dir,
                                  glue::glue("logrank_{group_df$subtype[i]}_EFS_cpgPLPstatus.RDS")))
  
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
  
  # take Ns from OS result so that Ns match
  efs_n <- sum(summary(km_os_plp$model)$table[,"records"])
  
  #add median EFS
  efs_stats[efs_stats$group == group_df$group[i],]$median_surv_years_plp <- summary(km_efs_plp$model)$table[,"median"][2]/365.25
  efs_stats[efs_stats$group == group_df$group[i],]$median_surv_years_no_plp <- summary(km_efs_plp$model)$table[,"median"][1]/365.25
  
  efs_diff_obj <- survdiff(survival::Surv(EFS_days, EFS_status) ~ cpgPLP_status,  
                          km_efs_plp$original_data)
  
  efs_stats[efs_stats$group == group_df$group[i],]$km_p <- 1 - pchisq(efs_diff_obj$chisq, length(efs_diff_obj$n) - 1)
  
  # extract EFS HRs and p-values
  efs_stats[efs_stats$group == group_df$group[i],]$HR <- exp(efs_df$estimate[efs_df$term == "cpgPLP_statusCPG P/LP"])
  efs_stats[efs_stats$group == group_df$group[i],]$p <- efs_df$p.value[efs_df$term == "cpgPLP_statusCPG P/LP"]
  
  efs_stats[efs_stats$group == group_df$group[i],]$CI_lower <- efs_ci_df["cpgPLP_statusCPG P/LP", "lower .95"]
  efs_stats[efs_stats$group == group_df$group[i],]$CI_upper <- efs_ci_df["cpgPLP_statusCPG P/LP", "upper .95"]
  
  # take Ns as the largest value in OS or EFS survival models 
  os_stats[os_stats$group == group_df$group[i],]$group_n <- max(c(os_n, efs_n))
  efs_stats[efs_stats$group == group_df$group[i],]$group_n <- max(c(os_n, efs_n))
  
}


# Some models have massively inflated p-values driven by single events; set these to NA
os_stats[abs(log10(os_stats$HR)) > 3 & !is.na(os_stats$HR), c("HR", "p", "CI_lower", "CI_upper")] <- NA
efs_stats[abs(log10(efs_stats$HR)) > 3 & !is.na(efs_stats$HR), c("HR", "p", "CI_lower", "CI_upper")] <- NA

# merge OS and EFS stats
survival_stats <- os_stats %>%
  bind_rows(efs_stats) 

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

survival_stats <- survival_stats %>%
  dplyr::mutate(group_plus_n = glue::glue("{group} (N={group_n})")) %>%
  dplyr::mutate(group_plus_n = fct_relevel(group_plus_n, 
                                           unique(group_plus_n)))


pdf(NULL)

# Save table to output
survival_stats %>% 
  dplyr::select(-p_label) %>%
write_tsv(file.path(results_dir, 
                    "median-survival-by-ancestry-cancer-group.tsv")) 

# plot median survival and KM p-values by histology/subtype
surv_type <- c("Overall survival", "Event-free survival")
names(surv_type) <- c("OS", "EFS")

for (i in 1:length(surv_type)){

  # generate KM p-value plot, ordering by p-values
  km_p_plot <- survival_stats %>%
    dplyr::filter(type == surv_type[i]) %>%
    arrange(desc(km_p)) %>%
    dplyr::mutate(group_plus_n = factor(group_plus_n,
                                             group_plus_n)) %>%
    
    ggplot(aes(x = -log10(km_p), y = factor(group_plus_n))) +
    geom_point(size = 3, color = "#00A087FF",
               show.legend = FALSE) + 
    labs(x = paste0("log10(", names(surv_type[i]), " K-M p) \n "), y = NULL) + 
    xlim(c(0, 3)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    theme_Publication() +
    theme(axis.text.x = ggtext::element_markdown())
  
  # plot median survival by P/LP carrier status and group, when median surv is reached
  median_surv_plot <- survival_stats %>%
    dplyr::filter(type == surv_type[i]) %>%
    pivot_longer(
      cols = c(median_surv_years_plp, median_surv_years_no_plp),
      names_to = "plp_status",
      names_prefix = "median_surv_years_",
      values_to = "median_surv_years"
    ) %>%
    dplyr::mutate(plp_status = case_when(
      plp_status == "plp" ~ "CPG P/LP carrier",
      TRUE ~ "non-carrier"
    )) %>%
    arrange(desc(km_p)) %>%
    dplyr::mutate(group_plus_n = factor(group_plus_n,
                                        unique(group_plus_n))) %>%
    
    ggplot(aes(x = median_surv_years, y = factor(group_plus_n),
               fill = plp_status)) +
    geom_col(col = "black",
           #  position = "dodge",
             position = position_dodge(width = 0.75),  # spacing within groups
             width = 0.7,
             show.legend = TRUE) + 
    labs(x = glue::glue("Median\n {surv_type[i]}"), y = NULL,
         fill = NULL) + 
    theme_Publication() +
    theme(axis.text.x = ggtext::element_markdown()) +
    theme(axis.text.y = element_blank(),
          legend.position = "top") +
    guides(fill = guide_legend(ncol = 2))
  
  # merge plots
  merged_plot <- ggarrange(km_p_plot,
                              median_surv_plot,
                              align = "h",
                              common.legend = TRUE,
                              widths = c(0.675, 0.3))
  
  # write to output
  ggsave(file.path(plot_dir, 
                   glue::glue("median-{names(surv_type[i])}-sig-byGroup.pdf")),
         merged_plot,
         width = 7, height = 7)
  
}

sessionInfo()
