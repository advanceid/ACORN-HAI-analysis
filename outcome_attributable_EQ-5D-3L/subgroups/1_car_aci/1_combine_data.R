# Clear 
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr, 
                 dplyr, 
                 labelled)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")

# Check for NA values in critical columns
columns_to_check <- c(
  "eq_5d_3l_1", "eq_5d_3l_2", "eq_5d_3l_3",
  "eq_5d_3l_4", "eq_5d_3l_5"
)
df_clean <- df[complete.cases(df[, columns_to_check]), ]

# Universal function to filter and clean data
process_data_subset <- function(df, filter_column) {
  if (!filter_column %in% colnames(df)) return(NULL)
  df <- df %>% filter(!is.na(.data[[filter_column]]))
  
  # Subject more than 12 years old
  df <- df %>% filter(age_new >= 12)
  
  # Relevel factors
  df$pathogen_combined_types <- droplevels(df$pathogen_combined_types)
  
  df$icu_hd_ap <- relevel(df$icu_hd_ap, ref = "No")
  
  
  # Set variable labels
  labels <- list(
    sex = "Sex",
    age_new = "Age (years)",
    country_region = "Region",
    country_income = "World Bank income status",
    hpd_admreason = "Primary admission reason",
    icu_hd_ap = "Admission to ICU/HD at enrollment",
    comorbidities_CCI = "Charlson comorbidity index",
    sofa_score = "SOFA score",
    severity_score_scale = "Severity score of disease",
    fbis_score = "FBIS score",
    pitt_score = "PITT score",
    qpitt_score = "qPITT score"
  )
  
  for (col in names(labels)) {
    if (col %in% colnames(df)) {
      df <- labelled::set_variable_labels(df, !!sym(col) := labels[[col]])
    }
  }
  
  return(df)
}

# Function to relevel factors ending with "-susceptible"
relevel_factors_by_suffix <- function(df) {
  factor_columns <- colnames(df)[sapply(df, is.factor)] 
  for (col in factor_columns) {
    levels_col <- levels(df[[col]])
    if (any(grepl("-susceptible$", levels_col))) { 
      ref_level <- grep("-susceptible$", levels_col, value = TRUE)
      if (length(ref_level) == 1) {
        df[[col]] <- relevel(df[[col]], ref = ref_level)
      }
    }
  }
  return(df)
}

# Filter and clean subsets for different pathogens
pathogen_columns <- c("aci_car", "ent_thir", "ent_car")
df_used <- lapply(pathogen_columns, function(col) process_data_subset(df_clean, col))

# Apply releveling to each subset in df_used
df_used <- lapply(df_used, function(sub_df) {
  if (!is.null(sub_df)) {
    sub_df <- relevel_factors_by_suffix(sub_df)
  }
  return(sub_df)
})

# Save data
saveRDS(df_used, "data/att_eq.RData")
###