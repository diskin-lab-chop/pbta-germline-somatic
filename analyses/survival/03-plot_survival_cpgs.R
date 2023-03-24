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
            "Neurofibroma plexiform",
            "CRANIO, ADAM", "EPN, PF A",
            "DMG, H3 K28", "HGG, H3 WT",
            "GNG/GNT, Other alteration", 
            "LGG, BRAF fusion", "LGG, Other alteration",
            "LGG, WT",
            "MB, Group3", "MB, Group 4", "MB, SHH"
)

file_names <- c("atrt", "cranio", "cpt", "epn", "hgg", "lgg", "mb", "mng", "mes", "gng_gnt", "nfp", 
                "cranio_ADAM", "epn_PF A", "hgg_k28", "hgg_WT", "gng_gnt_Other",
                "lgg_BRAF", "lgg_Other", "lgg_WT", "mb_Group3", "mb_Group4", "mb_SHH")
names(file_names) <- groups

# Define plot titles
titles <- c("ATRT", "Craniopharyngioma", "Choroid plexus tumor", 
            "Ependymoma", "HGG", "LGG", "Medulloblastoma", 
            "Meningioma", "Mesenchymal tumor",
            "Mixed neuronal-glial tumor", "Neurofibroma Plexiform", 
            "CRANIO, ADAM", "EPN, PF A",
            "DMG, H3 K28", "HGG, H3 WT",
            "GNG/GNT Other alteration", 
            "LGG, BRAF fusion", "LGG, Other alteration",
            "LGG, WT",
            "MB, Group3", "MB, Group 4", "MB_SHH")
names(titles) <- groups

# loop through cancer groups and subtypes and generate KM survival curves
for (group in groups){
  
  # read in km files
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_OS_cpgPLPstatus.RDS")
    )
  )
  
  km_efs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_EFS_cpgPLPstatus.RDS")
    )
  )
  
  # define output file
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{file_names[group]}_cpgPLPstatus.pdf"))
  
  # generate plot
  km_plot <- plotKM(model = list(km_os_result, km_efs_result),
                    variable = "cpgPLP_status",
                    combined = T, 
                    title = titles[group])
  
  # save plot
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 6, units = "in",
         device = "pdf")
  
  
}

# re-define `EPN, PF A`, as file name changes between km and coxph models
file_names["EPN, PF A"] <- "EPN_PF-A"

# loop through cancer groups and subtypes to generate OS forest plots
for (group in groups){
  
  # Read in coxph models
  if (grepl("Low-grade glioma|LGG", group)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  # define forest plot output file
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_OS_subtype_cpgPLPstatus.pdf"))
  
  # generate forest plot
  forest_plot <- plotForest(survival_result)
  
  # save plot
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}

# repeat plotting loop for EFS forest plots
for (group in groups){
  
  if (grepl("Low-grade glioma|LGG", group)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_resection_cpgPLPstatus.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_cpgPLPstatus.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_EFS_subtype_cpgPLPstatus.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}



