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

# Working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_eq.RData")
df_all <- as.data.frame(df_all[[1]])
df_used <- split(df_all, df_all$infection_types)

# Load fonts
loadfonts()

# Container for results (3 slots, default NULL)
result_save <- vector("list", 3)

system.time({
  for (i in 1:2) {
    # Only analyze i == 1, skip others
    if (i != 1) {
      result_save[[i]] <- NULL
      next
    }
    
    df <- df_used[[i]]
    
    # Data frame to store regression results
    results <- data.frame(variable = character(), 
                          coefficient = numeric(), 
                          lower_ci = numeric(), 
                          upper_ci = numeric(), 
                          p.value = character(), 
                          stringsAsFactors = FALSE)
    
    # Delete: "icu_hd_ap"
    variables <- c("age_new", "sex", 
                   "country_region",
                   "country_income",
                   "severity_score_scale",
                   "pathogen_combined_types",
                   "aci_car")
    
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
    
    # Manually structured rows (your custom table layout)
    rows <- list(
      c("Characteristics", NA, NA, NA, "Adjusted coefficient (95%CI)", "p.value"),
      c("Age", NA, NA, NA, NA, NA),
      c(NA, results[1, c(2:4, 6, 5)]),
      c("Sex", NA, NA, NA, NA, NA),
      c("Female", NA, NA, NA, "Ref", NA),
      c("Male", results[2, c(2:4, 6, 5)]),
      c("Region", NA, NA, NA, NA, NA),
      c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
      c("South-East Asian Region", results[3, c(2:4, 6, 5)]),
      c("Western Pacific Region", results[4, c(2:4, 6, 5)]),
      c("World Bank income status", NA, NA, NA, NA, NA),
      c("High income", NA, NA, NA, "Ref", NA),
      c("Upper middle income", results[5, c(2:4, 6, 5)]),
      c("Lower middle income", results[6, c(2:4, 6, 5)]),
      c("Primary admission reason", NA, NA, NA, NA, NA),
      c("Infectious disease", NA, NA, NA, NA, NA),
      c("Gastrointestinal disorder", NA, NA, NA, NA, NA),
      c("Pulmonary disease", NA, NA, NA, NA, NA),
      c("Others", NA, NA, NA, NA, NA),
      c("Charlson comorbidity index", NA, NA, NA, NA, NA),
      c(NA, NA, NA, NA, NA, NA),
      c("Severity score of disease", NA, NA, NA, NA, NA),
      c(NA, results[7, c(2:4, 6, 5)]),
      c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
      c("No", NA, NA, NA, NA, NA),
      c("Yes", NA, NA, NA, NA, NA),
      c("Type of pathogens", NA, NA, NA, NA, NA),
      c("Monomicrobial", NA, NA, NA, "Ref", NA),
      c("Polymicrobial", results[8, c(2:4, 6, 5)]),
      c("Acinetobacter spp.", NA, NA, NA, NA, NA),
      c("Carbapenem-susceptible", NA, NA, NA, "Ref", NA),
      c("Carbapenem-resistant", results[9, c(2:4, 6, 5)])
    )
    
    result <- do.call(rbind, rows) %>% as.data.frame()
    result[, 5:6] <- lapply(result[, 5:6], as.character)
    
    result_save[[i]] <- result[, c(1, 5:6)]
  }
})

# Save
saveRDS(result_save, "data/attfbis_value_multi_3.RData")
