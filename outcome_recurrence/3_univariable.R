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
df_used <- result_save <- list()

system.time({
  for (i in 1:3) {
    df <- df_data[[i]]
    
    # Convert outcome variable to numeric
    df$event <- as.numeric(df$event)
    
    # Define explanatory variables
    explanatory <- c("age_new", "sex", 
                     "country_region", "country_income",
                     "hpd_admreason", "comorbidities_CCI", 
                     "severity_score_scale",
                     "icu_hd_ap", 
                     "aci_car", "ent_thir", "ent_car")
    
    # Iterate over explanatory variables
    results_list <- list()
    
    for (var in explanatory) {
      # Create model matrix for the current variable
      covariates <- as.matrix(model.matrix(as.formula(paste("~", var)), data = df)[,-1])
      
      # Fit the model
      fmodel <- cmprsk::crr(
        ftime = df$time,
        fstatus = df$event,
        cov1 = covariates,
        failcode = 1,
        cencode = 0
      )
      
      # Get summary for the current variable
      summary_table <- custom_summary_crr(fmodel, var)
      
      # Add variable levels to the summary table, excluding the reference level
      if (is.factor(df[[var]]) || is.character(df[[var]])) {
        levels_var <- levels(df[[var]])
        if (length(levels_var) > (nrow(summary_table) + 1)) {
          summary_table$level <- levels_var[2:(nrow(summary_table) + 1)]
        } else {
          summary_table$level <- levels_var[-1]
        }
      } else {
        summary_table$level <- as.character(NA) 
      }
      
      # Reorder columns to place 'level' as the second column
      summary_table <- summary_table %>% select(variable, level, everything())
      
      # Store results in the list
      results_list[[var]] <- summary_table
    }
    
    # Combine all summary tables into a single data frame
    combined_results <- bind_rows(results_list, .id = "variable")
    
    
    #
    hazard_ratio_df <- transform(combined_results, 
                                 Hazard_Ratio = sprintf("%.3f", exp_coef), 
                                 Lower_CI = sprintf("%.3f", lower_95), 
                                 Upper_CI = sprintf("%.3f", upper_95))
    
    hazard_ratio_df$`Hazard ratio (95%CI)` = paste0(hazard_ratio_df$Hazard_Ratio, " [", hazard_ratio_df$Lower_CI, ", ", hazard_ratio_df$Upper_CI, "]", sep = "")
    
    df_used <- hazard_ratio_df %>% 
      select(variable, level, Hazard_Ratio, Lower_CI, Upper_CI, `Hazard ratio (95%CI)`, p_value)
    
    # Save HR
    rows <- list(
      c("Characteristics", NA, NA, NA, "Crude HR (95%CI)", "p.value"),
      c("Age", NA, NA, NA, NA, NA),
      c(df_used[1, c(2:7)]),
      c("Sex", NA, NA, NA, NA, NA),
      c("Female", NA, NA, NA, "Ref", NA),
      c(df_used[2, c(2:7)]),
      c("Region", NA, NA, NA, NA, NA),
      c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
      c(df_used[3, c(2:7)]),
      c(df_used[4, c(2:7)]),
      c("World Bank income status", NA, NA, NA, NA, NA),
      c("High income", NA, NA, NA, "Ref", NA),
      c(df_used[5, c(2:7)]),
      c(df_used[6, c(2:7)]),
      c("Primary admission reason", NA, NA, NA, NA, NA),
      c("Infectious disease", NA, NA, NA, "Ref", NA),
      c(df_used[7, c(2:7)]),
      c(df_used[8, c(2:7)]),
      c(df_used[9, c(2:7)]),
      c("Charlson comorbidity index", NA, NA, NA, NA, NA),
      c(df_used[10, c(2:7)]),
      c("Severity score of disease", NA, NA, NA, NA, NA),
      c(df_used[11, c(2:7)]),
      c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
      c("No", NA, NA, NA, "Ref", NA),
      c(df_used[12, c(2:7)]),
      c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c(df_used[13, c(2:7)]),
      c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c(df_used[14, c(2:7)]),
      c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c(df_used[15, c(2:7)])
    )
    
    # Bind the rows into a matrix
    result <- do.call(rbind, rows)
    
    #
    result <- as.data.frame(result)
    
    result[, 5:6] <- lapply(result[, 5:6], as.character)
    result[, 2:4] <- lapply(result[, 2:4], as.numeric)
    
    result_save[[i]] <- result[, c(1, 5, 6)]
   
  }
})

# Save
saveRDS(result_save, "data/recur_table_univariable.RData")
###