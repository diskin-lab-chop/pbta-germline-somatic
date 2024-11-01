# Functions for generating summary statistics 
#
# Ryan Corbett
#
# 2023

# Function to stratify patient cohort by supplied variable, and generate numnber of frequency of patients in each group

summarize_freq <- function(hist, var){
  
  freq <- hist %>%
    count(!!rlang::ensym(var)) %>%
    dplyr::mutate(`Percent Cohort` = round(n/nrow(hist)*100, 2)) %>%
    dplyr::rename("Number Cohort" = "n") %>%
    dplyr::mutate(`Number Cohort` = glue::glue("{`Number Cohort`} ({`Percent Cohort`}%)")) %>%
    dplyr::select(-`Percent Cohort`) %>%
    column_to_rownames(var) %>%
    t() %>%
    as.data.frame()
  
  return(freq)
  
}


# Function to stratify patient cohort by two supplied variables, and generate number of patients in each group

summarize_count <- function(hist, var1, var2){
  
  summary <- hist %>%
    count(!!rlang::ensym(var1), !!rlang::ensym(var2)) %>%
    spread(!!rlang::ensym(var2), n) 
  
  return(summary)
}


plp_enrichment <- function(hist, var){
  
  total_plp <- sum(hist$cpgPLP_status == "CPG P-LP")
  total_no_plp <- sum(hist$cpgPLP_status == "No CPG P-LP")
  
  enr_df <- hist %>%
    ## group by junction and calculate means
    select(!!rlang::ensym(var), cpgPLP_status) %>%
    group_by(!!rlang::ensym(var), cpgPLP_status) %>%
    count() %>%
    ungroup() %>%
    # Spread to wide format to get separate columns for "High" and "Low"
    pivot_wider(names_from = cpgPLP_status, values_from = n, values_fill = list(n = 0)) %>%
    rowwise() %>%
    mutate(
      Fisher_Test = list(
        fisher.test(
          matrix(
            c(`CPG P-LP`, total_plp - `CPG P-LP`,  # Counts of High and the absence of High
              `No CPG P-LP`, total_no_plp - `No CPG P-LP`),    # Counts of Low and the absence of Low
            nrow = 2
          )
        )
      ),
      OR = Fisher_Test$estimate,
      P_Value = Fisher_Test$p.value,
      CI_lower = Fisher_Test$conf.int[1],
      CI_upper = Fisher_Test$conf.int[2]
    ) %>%
    select(!!rlang::ensym(var), `CPG P-LP`, `No CPG P-LP`, OR, P_Value,
           CI_lower, CI_upper) %>%
    ungroup() %>%
    dplyr::mutate(adjusted_p = p.adjust(P_Value, method = "fdr")) %>%
    arrange(desc(OR)) %>%
    dplyr::mutate(fdr_label = case_when(
      adjusted_p < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    dplyr::mutate(OR = case_when(
      OR > 100 ~ 100, 
      OR < 0.01 ~ 0.01, 
      TRUE ~ OR
    )) %>%
    dplyr::mutate(CI_upper = case_when(
      is.infinite(CI_upper) ~ 10000, 
      TRUE ~ CI_upper
    )) %>%
    dplyr::mutate(CI_lower = case_when(
      CI_lower == 0 ~ 0.0001, 
      TRUE ~ CI_lower
    ))
  
  return(enr_df)
}

