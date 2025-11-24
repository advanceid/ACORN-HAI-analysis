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

# Load data
df_data_ori <- readRDS("data/clean_data.RData")

# Split BSI
df_data <- list(
  df_data_ori[[1]],
  split(df_data_ori[[2]], df_data_ori[[2]]$infection_types)[[2]],
  split(df_data_ori[[2]], df_data_ori[[2]]$infection_types)[[3]]
)

# ----------------------------------------------
# Define custom summary function for crr model
custom_summary_crr <- function(model, var) {
  # Extract coefficients and variance-covariance matrix
  coefs <- model$coef
  var_matrix <- model$var
  
  # Compute standard errors, z-values, and p-values
  se <- sqrt(diag(var_matrix))
  z_values <- coefs / se
  p_values <- 2 * pnorm(-abs(z_values))
  
  # Format p-values
  p_values_formatted <- ifelse(p_values < 0.001, "<0.001",
                               ifelse(p_values < 0.01, sprintf("%.3f", p_values), 
                                      sprintf("%.2f", p_values)))
  
  # Exponentiate coefficients and compute confidence intervals
  exp_coefs <- exp(coefs)
  conf_int <- matrix(NA, nrow = length(coefs), ncol = 2)
  for (i in 1:length(coefs)) {
    conf_int[i,] <- exp(c(coefs[i] - 1.96 * se[i], coefs[i] + 1.96 * se[i]))
  }
  
  # Create summary table
  summary_table <- data.frame(
    variable = var,
    coef = coefs,
    exp_coef = exp_coefs,
    se_coef = se,
    z = z_values,
    p_value = p_values_formatted,
    exp_neg_coef = 1 / exp_coefs,
    lower_95 = conf_int[,1],
    upper_95 = conf_int[,2]
  )
  
  return(summary_table)
}
# ----------------------------------------------
# Initialize list to store results
result_save <- list()

system.time({
  for (i in 1:3) {
    df <- df_data[[i]]
    
    if (i == 1) {
      result_save[[i]] <- NULL
      next
    }
    
    # Convert outcome variable to numeric
    df$event <- as.numeric(df$event)
    
    # x.factors
    if (i == 2) {
      # For factors, delete hpd_admreason, comorbidities_CCI
      x.factors <- model.matrix(~ df$sex + 
                                  df$country_region + df$country_income +
                                  df$icu_hd_ap +
                                  df$aci_car + df$ent_thir + df$ent_car)[,-1]
      
      # For num
      x.factors <- as.matrix(data.frame(x.factors, df$age_new,
                                        df$severity_score_scale))
    } else if (i == 3) {
      # For factors, delete comorbidities_CCI
      x.factors <- model.matrix(~ df$sex + 
                                  df$country_income +
                                  df$icu_hd_ap +
                                  df$aci_car + df$ent_thir + df$ent_car)[,-1]
      
      # For num
      x.factors <- as.matrix(data.frame(x.factors, df$age_new))
    }
    
    # Competing risk model
    fmodel <- cmprsk::crr(ftime = df$time, fstatus = df$event, x.factors,
                          failcode = 1, cencode = 0)
    
    # Summary model
    combined_results <- custom_summary_crr(fmodel, colnames(x.factors))
    
    #
    hazard_ratio_df <- transform(combined_results, 
                                 Hazard_Ratio = sprintf("%.3f", exp_coef), 
                                 Lower_CI = sprintf("%.3f", lower_95), 
                                 Upper_CI = sprintf("%.3f", upper_95))
    
    hazard_ratio_df$`Hazard ratio (95%CI)` = paste0(hazard_ratio_df$Hazard_Ratio, " [", hazard_ratio_df$Lower_CI, ", ", hazard_ratio_df$Upper_CI, "]", sep = "")
    
    df_used <- hazard_ratio_df %>% 
      select(variable, Hazard_Ratio, Lower_CI, Upper_CI, `Hazard ratio (95%CI)`, p_value)

  
    # Save HR
    if (i == 2) {
      rows <- list(
        c("Characteristics", NA, NA, NA, "Adjusted HR (95%CI)", "p.value"),
        c("Age", NA, NA, NA, NA, NA),
        c(NA, df_used[10, c(2:6)]),
        c("Sex", NA, NA, NA, NA, NA),
        c("Female", NA, NA, NA, "Ref", NA),
        c(NA, df_used[1, c(2:6)]),
        c("Region", NA, NA, NA, NA, NA),
        c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
        c("South-East Asian Region", df_used[2, c(2:6)]),
        c("Western Pacific Region", df_used[3, c(2:6)]),
        c("World Bank income status", NA, NA, NA, NA, NA),
        c("High income", NA, NA, NA, "Ref", NA),
        c("Upper middle income", df_used[4, c(2:6)]),
        c("Lower middle income", df_used[5, c(2:6)]),
        c("Primary admission reason", NA, NA, NA, NA, NA),
        c("Infectious disease", NA, NA, NA, NA, NA),
        c("Gastrointestinal disorder", NA, NA, NA, NA, NA),
        c("Pulmonary disease", NA, NA, NA, NA, NA),
        c("Others", NA, NA, NA, NA, NA),
        c("Charlson comorbidity index", NA, NA, NA, NA, NA),
        c(NA, NA, NA, NA, NA, NA),
        c("Severity score of disease", NA, NA, NA, NA, NA),
        c(NA, df_used[11, c(2:6)]),
        c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
        c("No", NA, NA, NA, "Ref", NA),
        c("Yes", df_used[6, c(2:6)]),
        c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
        c("Absent", NA, NA, NA, "Ref", NA),
        c("Present", df_used[7, c(2:6)]),
        c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
        c("Absent", NA, NA, NA, "Ref", NA),
        c("Present", df_used[8, c(2:6)]),
        c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
        c("Absent", NA, NA, NA, "Ref", NA),
        c("Present", df_used[9, c(2:6)])
      )
    } else if (i == 3) {
      rows <- list(
        c("Characteristics", NA, NA, NA, "Adjusted HR (95%CI)", "p.value"),
        c("Age", NA, NA, NA, NA, NA),
        c(NA, df_used[8, c(2:6)]),
        c("Sex", NA, NA, NA, NA, NA),
        c("Female", NA, NA, NA, "Ref", NA),
        c(NA, df_used[1, c(2:6)]),
        c("Region", NA, NA, NA, NA, NA),
        c("Eastern Mediterranean Region", NA, NA, NA, NA, NA),
        c("South-East Asian Region", NA, NA, NA, NA, NA),
        c("Western Pacific Region", NA, NA, NA, NA, NA),
        c("World Bank income status", NA, NA, NA, NA, NA),
        c("High income", NA, NA, NA, "Ref", NA),
        c("Upper middle income", df_used[2, c(2:6)]),
        c("Lower middle income", df_used[3, c(2:6)]),
        c("Primary admission reason", NA, NA, NA, NA, NA),
        c("Infectious disease", NA, NA, NA, NA, NA),
        c("Gastrointestinal disorder", NA, NA, NA, NA, NA),
        c("Pulmonary disease", NA, NA, NA, NA, NA),
        c("Others", NA, NA, NA, NA, NA),
        c("Charlson comorbidity index", NA, NA, NA, NA, NA),
        c(NA, NA, NA, NA, NA, NA),
        c("Severity score of disease", NA, NA, NA, NA, NA),
        c(NA, NA, NA, NA, NA, NA),
        c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
        c("No", NA, NA, NA, "Ref", NA),
        c("Yes", df_used[4, c(2:6)]),
        c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
        c("Absent", NA, NA, NA, "Ref", NA),
        c("Present", df_used[5, c(2:6)]),
        c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
        c("Absent", NA, NA, NA, "Ref", NA),
        c("Present", df_used[6, c(2:6)]),
        c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
        c("Absent", NA, NA, NA, "Ref", NA),
        c("Present", df_used[7, c(2:6)])
      )
    } 
  
    
    #
    result <- do.call(rbind, rows)
    result <- as.data.frame(result)
    
    result[, 5:6] <- lapply(result[, 5:6], as.character)
    result[, 2:4] <- lapply(result[, 2:4], as.numeric)
    
    # Save table
    result_save[[i]] <- result[, c(1, 5, 6)]
    
  }
})


# Save
saveRDS(result_save, "data/recur_table_multivariable_2.RData")
###
