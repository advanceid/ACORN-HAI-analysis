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

# Enterobacterales, carbapenem-resistant
df <- as.data.frame(df_all[[3]])

#
df$fbis_score <- factor(df$fbis_score, ordered = TRUE)

# Model, delete comorbidities_CCI
formula <- fbis_score ~
  age_new + sex + country_region + country_income +
  hpd_admreason + severity_score_scale + icu_hd_ap + 
  pathogen_combined_types +
  ent_car  

model <- MASS::polr(formula, data = df, method = "logistic", Hess = TRUE)

# Summarize model
model_summary <- summary(model)
p_model <- data.frame(model_summary$coefficients)
p_model$values <- p_model$Value
p_model$se_values <- p_model$Std..Error
p_model$t_values <- p_model$t.value

p_model$upper_ci <- p_model$values + 1.96 * p_model$se_values
p_model$lower_ci <- p_model$values - 1.96 * p_model$se_values

p_values_for_var <- 2 * (1 - pnorm(abs(p_model$t_values)))
p_model$p.value <- ifelse(p_values_for_var < 0.001, "<0.001",
                          ifelse(p_values_for_var < 0.01, sprintf("%.3f", p_values_for_var), 
                                 sprintf("%.2f", p_values_for_var)))

# Delete rows
del_row <- grepl("\\|", rownames(p_model))

p_model <- p_model[!del_row, ]

# Coef
coef_df <- transform(p_model, 
                     Coef = sprintf("%.2f", values), 
                     Lower_CI = sprintf("%.2f", lower_ci), 
                     Upper_CI = sprintf("%.2f", upper_ci))

coef_df$`coef (95%CI)` = paste0(coef_df$Coef, " [", coef_df$Lower_CI, ", ", coef_df$Upper_CI, "]", sep = "")

rows <- list(
  c("Characteristics", NA, NA, NA, "Adjusted coefficient (95%CI)", "p.value"),
  c("Age", NA, NA, NA, NA, NA),
  c(NA, coef_df[1, c(10:13, 9)]),
  c("Sex", NA, NA, NA, NA, NA),
  c("Female", NA, NA, NA, "Ref", NA),
  c("Male", coef_df[2, c(10:13, 9)]),
  c("Region", NA, NA, NA, NA, NA),
  c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
  c("South-East Asian Region", coef_df[3, c(10:13, 9)]),
  c("Western Pacific Region", coef_df[4, c(10:13, 9)]),
  c("World Bank income status", NA, NA, NA, NA, NA),
  c("High income", NA, NA, NA, "Ref", NA),
  c("Upper middle income", coef_df[5, c(10:13, 9)]),
  c("Lower middle income", coef_df[6, c(10:13, 9)]),
  c("Primary admission reason", NA, NA, NA, NA, NA),
  c("Infectious disease", NA, NA, NA, "Ref", NA),
  c("Gastrointestinal disorder", coef_df[7, c(10:13, 9)]),
  c("Pulmonary disease", coef_df[8, c(10:13, 9)]),
  c("Others", coef_df[9, c(10:13, 9)]),
  c("Charlson comorbidity index", NA, NA, NA, NA, NA),
  c(NA, NA, NA, NA, NA, NA),
  c("Severity score of disease", NA, NA, NA, NA, NA),
  c(NA, coef_df[10, c(10:13, 9)]),
  c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
  c("No", NA, NA, NA, "Ref", NA),
  c("Yes", coef_df[11, c(10:13, 9)]),
  c("Type of pathogens", NA, NA, NA, NA, NA),
  c("Monomicrobial", NA, NA, NA, "Ref", NA),
  c("Polymicrobial",coef_df[12, c(10:13, 9)]),
  c("Enterobacterales", NA, NA, NA, NA, NA),
  c("Carbapenem-susceptible", NA, NA, NA, "Ref", NA),
  c("Carbapenem-resistant", coef_df[13, c(10:13, 9)]))

# Bind the rows into a matrix
result <- do.call(rbind, rows) %>%
  as.data.frame()

result[, 5:6] <- lapply(result[, 5:6], as.character)
result[, 2:4] <- lapply(result[, 2:4], as.numeric)


result_save <- result[, c(1, 5, 6)]

# Save table
saveRDS(result_save, "data/attfbis_table_multi_1_overall.RData")
