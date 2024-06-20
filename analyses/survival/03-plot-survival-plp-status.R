# R. Corbett 2023
#
# Makes a pdf panels for PBTA-germline Kaplan-Meier survival analyses 

library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival")

# call survival functions 
source(file.path(analysis_dir, "util", "survival_models.R"))

# Declare input and plot directory
input_dir <- file.path(root_dir, "analyses", "survival", "results")
plots_dir <- file.path(root_dir, "analyses", "survival", "plots")
results_dir <- file.path(root_dir, "analyses", "survival", "results")

if(!interactive()) pdf(NULL)

subtype_file <- file.path(results_dir, "subtypes-for-survival.tsv")

# Define cancer groups and subtypes for plotting
groups <- c("Atypical Teratoid Rhabdoid Tumor",
            "Craniopharyngioma",
            "Choroid plexus tumor",
            "DIPG or DMG|Other high-grade glioma",
            "Ependymoma",
            "Low-grade glioma",
            "Medulloblastoma",
            "Meningioma",
            "Mesenchymal tumor",
            "Mixed neuronal-glial tumor",
            "Neurofibroma plexiform")


dir_names <- c("ATRT", "CRANIO", "CPT", "HGG", "EPN", "LGG", "MB", "MNG", "MES", "GNG-GNT", "NFP")
names(dir_names) <- groups

# Define plot titles
titles <- c("ATRT", "Craniopharyngioma", "Choroid plexus tumor", 
            "HGG", "Ependymoma", "LGG", "Medulloblastoma", 
            "Meningioma", "Mesenchymal tumor",
            "Mixed neuronal-glial tumor", "Neurofibroma Plexiform")
names(titles) <- groups

# loop through cancer groups and subtypes and generate KM survival curves
for (group in groups){
  
  input_dir <- file.path(analysis_dir, "results", dir_names[group])
  plots_dir <- file.path(analysis_dir, "plots", dir_names[group])
  
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
    
  }
  
  # read in km files
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{dir_names[group]}_OS_cpgPLPstatus.RDS")
    )
  )
  
  km_efs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{dir_names[group]}_EFS_cpgPLPstatus.RDS")
    )
  )
  
  # define output file
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{dir_names[group]}_cpgPLPstatus.pdf"))
  
  pdf(NULL)
  
  # generate plot
  km_plot <- plotKM(model = list(km_os_result, km_efs_result),
                    variable = "cpgPLP_status",
                    combined = T, 
                    title = titles[group])
  
  # save plot
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 6, units = "in",
         device = "pdf")
  
  
  # Read in coxph models
  if (grepl("Low-grade glioma|LGG", group)){
    os_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{dir_names[group]}_OS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    os_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{dir_names[group]}_OS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  if (sum(is.na(broom::tidy(os_result)$estimate)) != nrow(broom::tidy(os_result))) {
  
    # define forest plot output file
    os_forest_pdf <- file.path(plots_dir, 
                            glue::glue("forest_{dir_names[group]}_OS_subtype_cpgPLPstatus.pdf"))
    
    pdf(NULL)
    
    # generate forest plot
    os_forest_plot <- plotForest(os_result)
    
    # save plot
    ggsave(os_forest_pdf, os_forest_plot, width = 8, height = 3)
    
  }
  
  if (grepl("Low-grade glioma|LGG", group)){
    efs_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{dir_names[group]}_EFS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    efs_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{dir_names[group]}_EFS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  if (sum(is.na(broom::tidy(efs_result)$estimate)) != nrow(broom::tidy(efs_result))) {
    
    efs_forest_pdf <- file.path(plots_dir, 
                            glue::glue("forest_{dir_names[group]}_EFS_subtype_cpgPLPstatus.pdf"))
    
    
    efs_forest_plot <- plotForest(efs_result)
    
    ggsave(efs_forest_pdf, efs_forest_plot, width = 8, height = 3)
    
  }
  
}


subtype_df <- read_tsv(subtype_file)


for (i in 1:nrow(subtype_df)){
  
  input_dir <- file.path(analysis_dir, "results", subtype_df$hist[i])
  plots_dir <- file.path(analysis_dir, "plots", subtype_df$hist[i])
  
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
    
  }
  
  # read in km files
  subtype_km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_OS_cpgPLPstatus.RDS")
    )
  )
  
  subtype_km_efs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_EFS_cpgPLPstatus.RDS")
    )
  )
  
  # define output file
  subtype_km_output_pdf <- file.path(plots_dir, glue::glue("km_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_cpgPLPstatus.pdf"))
  
  pdf(NULL)
  
  # generate plot
  subtype_km_plot <- plotKM(model = list(subtype_km_os_result, subtype_km_efs_result),
                    variable = "cpgPLP_status",
                    combined = T, 
                    title = subtype_df$mol_sub_group[i])
  
  # save plot
  ggsave(subtype_km_output_pdf, subtype_km_plot,
         width = 10, height = 6, units = "in",
         device = "pdf")
  
  
  # Read in coxph models
  if (grepl("Low-grade glioma", subtype_df$plot_group[i])){
    subtype_os_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_OS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    subtype_os_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_OS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  if (sum(is.na(broom::tidy(subtype_os_result)$estimate)) != nrow(broom::tidy(subtype_os_result))) {
  
  # define forest plot output file
  subtype_os_forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_OS_subtype_cpgPLPstatus.pdf"))
  
  pdf(NULL)
  
  # generate forest plot
  subtype_os_forest_plot <- plotForest(subtype_os_result)
  
  # save plot
  ggsave(subtype_os_forest_pdf, subtype_os_forest_plot, width = 8, height = 3)
  
  }
  
  if (grepl("Low-grade glioma", subtype_df$plot_group[i])){
    subtype_efs_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_EFS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    subtype_efs_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_EFS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  if (sum(is.na(broom::tidy(subtype_efs_result)$estimate)) != nrow(broom::tidy(subtype_efs_result))) {
  
    subtype_efs_forest_pdf <- file.path(plots_dir, 
                            glue::glue("forest_{subtype_df$hist[i]}_{subtype_df$subtype[i]}_EFS_subtype_cpgPLPstatus.pdf"))
    
    
    subtype_efs_forest_plot <- plotForest(subtype_efs_result)
    
    ggsave(subtype_efs_forest_pdf, subtype_efs_forest_plot, width = 8, height = 3)
    
    dev.off()
    
  }
  
}
