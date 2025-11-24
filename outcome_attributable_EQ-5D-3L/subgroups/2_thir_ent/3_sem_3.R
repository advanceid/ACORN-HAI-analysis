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

# Define directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_eq.RData")

df_all <- as.data.frame(df_all[[2]])

# Split by infection types
df_used <- split(df_all, df_all$infection_types)

#
result_save <- S_coef_values_all <- factor_load_eq_all <- list()

for (i in 1:3){
  df <- df_used[[i]]
  
  # Ensure `ent_thir` is numeric 
  df$ent_thir <- ifelse(grepl("resistant$", df$ent_thir), 1, 0)
  
  if (i == 1 || i == 2){
    result_save[[i]] <- NULL
    next
  }
  
  # Define SEM model
 if (i == 3) {
   sem_model <- '
   # 1. Define latent variable eq_5d_3l_latent
   eq_5d_3l_latent =~ eq_5d_3l_1 + eq_5d_3l_4 + eq_5d_3l_5
   out =~ eq_5d_3l_2 + eq_5d_3l_3

   # 2. Structural equation model, delete comorbidities_CCI
   eq_5d_3l_latent ~ age_new + sex + ent_thir
   out ~ eq_5d_3l_latent
  '
 }
  
  # Run SEM
  fit <- sem(sem_model, data = df, se = "bootstrap", bootstrap = 1000)
  
  # CI
  ci_boot <- standardizedSolution_boot_ci(fit)
  param_est <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  
  # Coef
  important_vars <- c("age_new", "sex", "comorbidities_CCI", "severity_score_scale", "icu_hd_ap", "pathogen_combined_types", "ent_thir")
  
  S_coef_values <- ci_boot %>%
    filter(lhs == "eq_5d_3l_latent" & rhs %in% important_vars)%>%
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
  
  # Factor load
  factor_load_eq <- ci_boot %>%
    filter(lhs %in% c("eq_5d_3l_latent", "out") & grepl("^eq_5d_3l_[1-5]$", rhs))%>%
    select(rhs, est.std, boot.ci.lower, boot.ci.upper, pvalue) %>%
    mutate(pvalue = ifelse(pvalue < 0.001, "<0.001",
                           ifelse(pvalue < 0.01, sprintf("%.3f", pvalue), 
                                  sprintf("%.2f", pvalue)))) %>%
    mutate(
      est.std = round(est.std, 3),
      boot.ci.lower = round(boot.ci.lower, 3),
      boot.ci.upper = round(boot.ci.upper, 3)
    )
  
  #
  if (i == 3) {
    rows <- list(
      c("Characteristics", NA, NA, NA, "Unstandardized coefficient (95%CI)", "p.value", NA, NA, NA, "Standardized coefficient (95%CI)"),
      c("Age", rep(NA,9)),
      c(NA, coef_df[1, c(2:4, 6, 5, 7:10)]),
      c("Sex", rep(NA,9)),
      c("Female", NA, NA, NA, "Ref", rep(NA,5)),
      c("Male", coef_df[2, c(2:4, 6, 5, 7:10)]),
      c("Charlson comorbidity index", rep(NA,9)),
      c(NA, rep(NA,9)),
      c("Severity score of disease", rep(NA,9)),
      c(NA, rep(NA,9)),
      c("Admission to ICU/HD at enrollment", rep(NA,9)),
      c("No", rep(NA,9)),
      c("Yes", rep(NA,9)),
      c("Type of pathogens", rep(NA,9)),
      c("Monomicrobial", rep(NA,9)),
      c("Polymicrobial", rep(NA,9)),
      c("Enterobacterales", rep(NA,9)),
      c("Third-generation cephalosporin-susceptible", NA, NA, NA, "Ref", rep(NA,5)),
      c("Third-generation cephalosporin-resistant", coef_df[3, c(2:4, 6, 5, 7:10)]))
  }
  
  # Bind the rows into a matrix
  result <- do.call(rbind, rows)
  
  #
  result <- as.data.frame(result)
  
  result[, c(5:6, 10)] <- lapply(result[, c(5:6, 10)], as.character)
  result[, c(2:4, 7:9)] <- lapply(result[, c(2:4, 7:9)], as.numeric)
  
  result_save[[i]] <- result
  
  #
  # Rename
  S_coef_values <- S_coef_values %>%
    mutate(rhs = case_when(
      rhs == "age_new" ~ "Age (years)",
      rhs == "sex" ~ "Male (vs female)",
      rhs == "severity_score_scale" ~ "Severity score of disease",
      rhs == "icu_hd_ap" ~ "Admission to ICU/HD at enrollment (yes vs no)",
      rhs == "pathogen_combined_types" ~ "Polymicrobial (vs monomicrobial)",
      rhs == "ent_thir" ~ "3GCRE (vs 3GCSE)",
      TRUE ~ rhs
    ))
  
  factor_load_eq <- factor_load_eq %>%
    mutate(rhs = case_when(
      rhs == "eq_5d_3l_1" ~ "Mobility",
      rhs == "eq_5d_3l_2" ~ "Self-care",
      rhs == "eq_5d_3l_3" ~ "Usual activities",
      rhs == "eq_5d_3l_4" ~ "Pain/discomfort",
      rhs == "eq_5d_3l_5" ~ "Anxiety/depression",
      TRUE ~ rhs
    ))
  
  S_coef_values_all[[i]] <- S_coef_values
  factor_load_eq_all[[i]] <- factor_load_eq
}

# Save
saveRDS(result_save, "data/attfbis_table_multi_3_sub.RData")
saveRDS(factor_load_eq_all, "data/attfbis_factor_load_3_sub.RData")
saveRDS(S_coef_values_all, "data/attfbis_S_coef_3_sub.RData")