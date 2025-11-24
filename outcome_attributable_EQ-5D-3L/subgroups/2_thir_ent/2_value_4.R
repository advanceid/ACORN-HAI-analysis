# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 glmnet,
                 car,
                 gt,
                 extrafont)
})
# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_eq.RData")

df_all <- as.data.frame(df_all[[2]])
df_used <- split(df_all, df_all$infection_types)

# Load fonts
loadfonts()

#
result_save <- list()

system.time({
  for (i in 1:3) {
    df <- df_used[[i]]
    
    if (i == 2){
      result_save[[i]] <- NULL
      next
    }
    
    # Initialize the data frame to store results
    results <- data.frame(variable = character(), 
                          coefficient = numeric(), 
                          lower_ci = numeric(), 
                          upper_ci = numeric(), 
                          p.value = character(), 
                          stringsAsFactors = FALSE)
    
    # Potential variables
    if (i == 1) {
      # delete "hpd_admreason"
      variables <- c("age_new", "sex", 
                     "country_income",
                     "comorbidities_CCI", 
                     "severity_score_scale",
                     "icu_hd_ap", "pathogen_combined_types",
                     "ent_thir")
    } else if (i == 3) {
      # delete "icu_hd_ap"
      variables <- c("age_new", "sex", 
                     "country_income",
                     "comorbidities_CCI", 
                     "severity_score_scale",
                     "pathogen_combined_types",
                     "ent_thir")
    } 
    
    # Model
    formula <- as.formula(paste("eq_5d_3l ~", paste(variables, collapse = " + ")))
    model <- glm(formula, family = gaussian(link = "identity"), data = df)
    sum_model <- as.data.frame(summary(model)$coefficients)
    
    # Skip intercept
    for (row in rownames(sum_model)) {
      if (row != "(Intercept)") {
        coef <- sum_model[row, "Estimate"]
        se <- sum_model[row, "Std. Error"]
        p <- sum_model[row, "Pr(>|t|)"]
        lower_ci <- coef - 1.96 * se
        upper_ci <- coef + 1.96 * se
        significance <- ifelse(p < 0.001, "<0.001", 
                               ifelse(p < 0.01, "<0.01", 
                                      ifelse(p < 0.05, "<0.05", sprintf("%.3f", p))))
        results <- rbind(results, data.frame(variable = row, 
                                             coefficient = coef, 
                                             lower_ci = lower_ci, 
                                             upper_ci = upper_ci, 
                                             p.value = significance))
      }
    }
    
    # Final data frame
    results <- results %>%
      mutate(across(c(coefficient, lower_ci, upper_ci), ~ sprintf("%.3f", .))) %>%
      mutate(`Coefficient (95%CI)` = paste0(coefficient, " [", lower_ci, ", ", upper_ci, "]"))
    
    # Rows
    if (i == 1) {
      rows <- list(
        c("Characteristics", NA, NA, NA, "Adjusted coefficient (95%CI)", "p.value"),
        c("Age", NA, NA, NA, NA, NA),
        c(NA, results[1, c(2:4, 6, 5)]),
        c("Sex", NA, NA, NA, NA, NA),
        c("Female", NA, NA, NA, "Ref", NA),
        c("Male", results[2, c(2:4, 6, 5)]),
        c("World Bank income status", NA, NA, NA, NA, NA),
        c("High income", NA, NA, NA, "Ref", NA),
        c("Upper middle income", results[3, c(2:4, 6, 5)]),
        c("Lower middle income", results[4, c(2:4, 6, 5)]),
        c("Primary admission reason", NA, NA, NA, NA, NA),
        c("Infectious disease", NA, NA, NA, NA, NA),
        c("Gastrointestinal disorder", NA, NA, NA, NA, NA),
        c("Pulmonary disease", NA, NA, NA, NA, NA),
        c("Others", NA, NA, NA, NA, NA),
        c("Charlson comorbidity index", NA, NA, NA, NA, NA),
        c(NA, results[5, c(2:4, 6, 5)]),
        c("Severity score of disease", NA, NA, NA, NA, NA),
        c(NA, results[6, c(2:4, 6, 5)]),
        c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
        c("No", NA, NA, NA, "Ref", NA),
        c("Yes", results[7, c(2:4, 6, 5)]),
        c("Type of pathogens", NA, NA, NA, NA, NA),
        c("Monomicrobial", NA, NA, NA, "Ref", NA),
        c("Polymicrobial", results[8, c(2:4, 6, 5)]),
        c("Enterobacterales", NA, NA, NA, NA, NA),
        c("Third-generation cephalosporin-susceptible", NA, NA, NA, "Ref", NA),
        c("Third-generation cephalosporin-resistant", results[9, c(2:4, 6, 5)])
      )
    } else if (i == 3) {
      rows <- list(
        c("Characteristics", NA, NA, NA, "Adjusted coefficient (95%CI)", "p.value"),
        c("Age", NA, NA, NA, NA, NA),
        c(NA, results[1, c(2:4, 6, 5)]),
        c("Sex", NA, NA, NA, NA, NA),
        c("Female", NA, NA, NA, "Ref", NA),
        c("Male", results[2, c(2:4, 6, 5)]),
        c("World Bank income status", NA, NA, NA, NA, NA),
        c("High income", NA, NA, NA, "Ref", NA),
        c("Upper middle income", results[3, c(2:4, 6, 5)]),
        c("Lower middle income", results[4, c(2:4, 6, 5)]),
        c("Primary admission reason", NA, NA, NA, NA, NA),
        c("Infectious disease", NA, NA, NA, NA, NA),
        c("Gastrointestinal disorder", NA, NA, NA, NA, NA),
        c("Pulmonary disease", NA, NA, NA, NA, NA),
        c("Others", NA, NA, NA, NA, NA),
        c("Charlson comorbidity index", NA, NA, NA, NA, NA),
        c(NA, results[5, c(2:4, 6, 5)]),
        c("Severity score of disease", NA, NA, NA, NA, NA),
        c(NA, results[6, c(2:4, 6, 5)]),
        c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
        c("No", NA, NA, NA, NA, NA),
        c("Yes", NA, NA, NA, NA, NA),
        c("Type of pathogens", NA, NA, NA, NA, NA),
        c("Monomicrobial", NA, NA, NA, "Ref", NA),
        c("Polymicrobial", results[7, c(2:4, 6, 5)]),
        c("Enterobacterales", NA, NA, NA, NA, NA),
        c("Third-generation cephalosporin-susceptible", NA, NA, NA, "Ref", NA),
        c("Third-generation cephalosporin-resistant", results[8, c(2:4, 6, 5)])
      )
    } 
    
    result <- do.call(rbind, rows) %>% as.data.frame()
    result[, 5:6] <- lapply(result[, 5:6], as.character)
    
    result_save[[i]] <- result[, c(1, 5:6)]
    
  }
})

# Save 
saveRDS(result_save, "data/attfbis_value_multi_2.RData")
