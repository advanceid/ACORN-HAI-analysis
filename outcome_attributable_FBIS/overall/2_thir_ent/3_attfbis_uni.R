# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 MASS)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_fbis.RData")

# Enterobacterales, third-generation cephalosporin-resistant
df <- as.data.frame(df_all[[2]])

# Univariable analysis
df$fbis_score <- factor(df$fbis_score, ordered = TRUE)

# potential variables
variables <- c("age_new", "sex", 
               "country_region", "country_income",
               "hpd_admreason", 
               "comorbidities_CCI", 
               "severity_score_scale",
               "icu_hd_ap", "pathogen_combined_types",
               "ent_thir")

# Initialize the data frame to store results
p_values <- data.frame(p.value = character(), ratio = numeric(), 
                       lower_ci = numeric(), upper_ci = numeric(), 
                       stringsAsFactors = FALSE)

# 
for (var in variables){
  formula <- as.formula(paste("fbis_score ~", var))
  univariate_model <- MASS::polr(formula, data = df, method = "logistic", Hess = TRUE)
  # Summary model
  model_summary <- summary(univariate_model)
  
  variable_names <- rownames(model_summary$coefficients)
  values <- model_summary$coefficients[, 1]
  se_values <- model_summary$coefficients[, 2]
  t_values <- model_summary$coefficients[, 3]
  
  p_values_for_var <- 2 * (1 - pnorm(abs(t_values)))
  
  for (i in seq_along(variable_names)) {
    if (!grepl("\\|", variable_names[i])) {
      p.values <- p_values_for_var[i]
      value <- values[i]
      se <- se_values[i]
      
      coef <- value
      upper_ci <- value + 1.96 * se
      lower_ci <- value - 1.96 * se
      
      significance <- ifelse(p.values < 0.001, "<0.001",
                             ifelse(p.values < 0.01, sprintf("%.3f", p.values), 
                                    sprintf("%.2f", p.values)))
      
      p_values <- rbind(p_values, data.frame(
        p.value = significance, 
        coef = coef, 
        lower_ci = lower_ci, 
        upper_ci = upper_ci
      ))
    }
  }
}
  
# Display the final data frame
p_values
  
# -------------------------
coef_df <- transform(p_values, 
                     Coef = sprintf("%.2f", coef), 
                     Lower_CI = sprintf("%.2f", lower_ci), 
                     Upper_CI = sprintf("%.2f", upper_ci))
  
coef_df$`coef (95%CI)` = paste0(coef_df$Coef, " [", coef_df$Lower_CI, ", ", coef_df$Upper_CI, "]", sep = "")
  
###
rows <- list(
    c("Characteristics", NA, NA, NA, "Crude coefficient (95%CI)", "p.value"),
    c("Age", NA, NA, NA, NA, NA),
    c(NA, coef_df[1, c(5:8, 1)]),
    c("Sex", NA, NA, NA, NA, NA),
    c("Female", NA, NA, NA, "Ref", NA),
    c("Male", coef_df[2, c(5:8, 1)]),
    c("Region", NA, NA, NA, NA, NA),
    c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
    c("South-East Asian Region", coef_df[3, c(5:8, 1)]),
    c("Western Pacific Region", coef_df[4, c(5:8, 1)]),
    c("World Bank income status", NA, NA, NA, NA, NA),
    c("High income", NA, NA, NA, "Ref", NA),
    c("Upper middle income", coef_df[5, c(5:8, 1)]),
    c("Lower middle income", coef_df[6, c(5:8, 1)]),
    c("Primary admission reason", NA, NA, NA, NA, NA),
    c("Infectious disease", NA, NA, NA, "Ref", NA),
    c("Gastrointestinal disorder", coef_df[7, c(5:8, 1)]),
    c("Pulmonary disease", coef_df[8, c(5:8, 1)]),
    c("Others", coef_df[9, c(5:8, 1)]),
    c("Charlson comorbidity index", NA, NA, NA, NA, NA),
    c(NA, coef_df[10, c(5:8, 1)]),
    c("Severity score of disease", NA, NA, NA, NA, NA),
    c(NA, coef_df[11, c(5:8, 1)]),
    c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
    c("No", NA, NA, NA, "Ref", NA),
    c("Yes", coef_df[12, c(5:8, 1)]),
    c("Type of pathogens", NA, NA, NA, NA, NA),
    c("Monomicrobial", NA, NA, NA, "Ref", NA),
    c("Polymicrobial", coef_df[13, c(5:8, 1)]),
    c("Enterobacterales", NA, NA, NA, NA, NA),
    c("Third-generation cephalosporin-susceptible", NA, NA, NA, "Ref", NA),
    c("Third-generation cephalosporin-resistant", coef_df[14, c(5:8, 1)])
  )
  
# Bind the rows into a matrix
result <- do.call(rbind, rows)
  
#
result <- as.data.frame(result)
  
result[, 5:6] <- lapply(result[, 5:6], as.character)
result[, 2:4] <- lapply(result[, 2:4], as.numeric)

result_save <- result[, c(1, 5, 6)]
  
# Save
saveRDS(result_save, "data/attfbis_table_univariable_overall.RData")
###