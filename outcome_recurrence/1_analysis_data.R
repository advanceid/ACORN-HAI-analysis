# Clear 
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, 
                 magrittr,
                 crrSC, 
                 cmprsk)
})

# Load the data
df_data <- readRDS("data/rp_data.RData")

# Initialize list to store results
df_used <- list()

system.time({
  for (i in 1:3) {
    df <- df_data[[i]]
    
    # Convert outcome variable to numeric
    df$event <- as.numeric(df$event)
    
    # Create a function to check if a variable's value ends with "-resistant"
    is_resistant <- function(x) {
      ifelse(grepl("-resistant$", x), 1, 0)
    }
    
    # CRA, 3GCRE, CRE, CRP, VRE, MRSA
    # Set factors
    cols <- c("aci_car", "ent_car", "ent_thir", "pse_car", "entc_van", "sa_meth")
    df[cols] <- lapply(df[cols], function(x) factor(is_resistant(x), levels = c(0, 1), labels = c("Absent", "Present")))
    
    # Define explanatory variables
    explanatory <- c("age_new", "sex", 
                     "country_region", "country_income",
                     "hpd_admreason", "comorbidities_CCI", 
                     "severity_score_scale",
                     "icu_hd_ap", "infection_types",
                     "aci_car", "ent_thir", "ent_car", 
                     "pse_car", "entc_van", "sa_meth")
    
    # Select relevant columns and remove rows with missing values
    if (i == 3) {
      df <- df %>%
        select(time, time_new, event, event_new, all_of(explanatory)) %>%
        filter(complete.cases(.), time > 0)
    } else {
      df <- df %>%
        select(time, event, all_of(explanatory)) %>%
        filter(complete.cases(.), time > 0)
    }
    
    # Set reference
    df$icu_hd_ap <- relevel(factor(df$icu_hd_ap), ref = "No")
    
    df_used[[i]] <- df
    
  }
})

# Save
saveRDS(df_used, "data/clean_data.RData")
###
