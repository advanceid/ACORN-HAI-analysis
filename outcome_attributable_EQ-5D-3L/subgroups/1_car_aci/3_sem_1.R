# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 fastDummies,
                 semhelpinghands,
                 lavaan)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_eq.RData")

df_all <- as.data.frame(df_all[[1]])

# Split by infection types
df_used <- split(df_all, df_all$infection_types)

#
result_save <- list()

for (i in 1:2){
  df <- df_used[[i]]
  
  # Ensure `aci_car` is numeric 
  df$aci_car <- ifelse(grepl("resistant$", df$aci_car), 1, 0)
  
  sem_model <- '
  # 1. Define latent variable eq_5d_3l_latent
  eq_5d_3l_latent =~ eq_5d_3l_1 + eq_5d_3l_4 + eq_5d_3l_5
  
  out =~ eq_5d_3l_2 + eq_5d_3l_3

  # 2. Structural equation model
  eq_5d_3l_latent ~ age_new + sex + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + aci_car
  out ~ eq_5d_3l_latent
  '
  
  # Run SEM
  fit <- sem(sem_model, data = df, se = "bootstrap", bootstrap = 1000)
  
  # CI
  ci_boot <- standardizedSolution_boot_ci(fit)
  param_est <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  
  # Coef
  important_vars <- c("age_new", "sex", "comorbidities_CCI", "severity_score_scale", "icu_hd_ap", "pathogen_combined_types", "aci_car")
  
  S_coef_values <- ci_boot %>%
    filter(lhs %in% c("eq_5d_3l_latent", "out") & rhs %in% important_vars)%>%
    select(rhs, est.std, boot.ci.lower, boot.ci.upper, pvalue) %>%
    mutate(pvalue = ifelse(pvalue < 0.001, "<0.001",
                           ifelse(pvalue < 0.01, sprintf("%.3f", pvalue), 
                                  sprintf("%.2f", pvalue)))) %>%
    mutate(
      est.std = round(est.std, 2),
      boot.ci.lower = round(boot.ci.lower, 2),
      boot.ci.upper = round(boot.ci.upper, 2)
    ) %>%
    mutate(`S_coef (95%CI)` = paste0(est.std, " [", boot.ci.lower, ", ", boot.ci.upper, "]"))
  
  U_coef_values <- param_est %>%
    filter(lhs == "eq_5d_3l_latent" & rhs %in% important_vars) %>%
    select(rhs, est, ci.lower, ci.upper, pvalue) %>%
    mutate(pvalue = ifelse(pvalue < 0.001, "<0.001",
                           ifelse(pvalue < 0.01, sprintf("%.3f", pvalue), 
                                  sprintf("%.2f", pvalue))))%>%
    mutate(
      est = round(est, 2),
      ci.lower = round(ci.lower, 2),
      ci.upper = round(ci.upper, 2)
    ) %>%
    mutate(`coef (95%CI)` = paste0(est, " [", ci.lower, ", ", ci.upper, "]"))
  
  coef_df <- left_join(U_coef_values, S_coef_values[,-5], by = "rhs")
  
  rows <- list(
    c("Characteristics", NA, NA, NA, "Unstandardized coefficient (95%CI)", "p.value", NA, NA, NA, "Standardized coefficient (95%CI)"),
    c("Age", rep(NA,9)),
    c(NA, coef_df[1, c(2:4, 6, 5, 7:10)]),
    c("Sex", rep(NA,9)),
    c("Female", NA, NA, NA, "Ref", rep(NA,5)),
    c("Male", coef_df[2, c(2:4, 6, 5, 7:10)]),
    c("Charlson comorbidity index", rep(NA,9)),
    c(NA, coef_df[3, c(2:4, 6, 5, 7:10)]),
    c("Severity score of disease", rep(NA,9)),
    c(NA, coef_df[4, c(2:4, 6, 5, 7:10)]),
    c("Admission to ICU/HD at enrollment", rep(NA,9)),
    c("No", NA, NA, NA, "Ref", rep(NA,5)),
    c("Yes", coef_df[5, c(2:4, 6, 5, 7:10)]),
    c("Type of pathogens", rep(NA,9)),
    c("Monomicrobial",  NA, NA, NA, "Ref", rep(NA,5)),
    c("Polymicrobial", coef_df[6, c(2:4, 6, 5, 7:10)]),
    c("Acinetobacter spp.", rep(NA,9)),
    c("Carbapenem-susceptible", NA, NA, NA, "Ref", rep(NA,5)),
    c("Carbapenem-resistant", coef_df[7, c(2:4, 6, 5, 7:10)]))
  
  # Bind the rows into a matrix
  result <- do.call(rbind, rows)
  
  #
  result <- as.data.frame(result)
  
  result[, c(5:6, 10)] <- lapply(result[, c(5:6, 10)], as.character)
  result[, c(2:4, 7:9)] <- lapply(result[, c(2:4, 7:9)], as.numeric)
  
  result_save[[i]] <- result
}

# Save 
saveRDS(result_save, "data/attfbis_table_multi_1_sub.RData")
