# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 survminer,
                 rstpm2,
                 scales)
})

# Set working directory
wd <- "./"
setwd(wd)

# Load data
pop <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
df_list <- readRDS("data/att_first28.RData")

# Define pathogen variables
pathogens <- c("aci_car", "ent_thir", "ent_car", "pse_car",
               "entc_van", "sa_meth", "mdr_gnb", "mdr", "amr")

# Define the tvc settings 
tvc_list <- list(
  # For pathogens[1] == "aci_car"
  list(
    age_new = 3,
    sex = 3,
    country_region = 3,
    pathogen_combined_types = 3,
    infection_types = 3
  ),
  # For pathogens[2] == "ent_thir"
  list(
    age_new = 3,
    country_region = 3,
    country_income = 3,
    severity_score_scale = 3,
    infection_types = 3
  ),
  # For pathogens[3] == "ent_car"
  list(
    age_new = 3,
    country_income = 3,
    hpd_admreason = 3,
    severity_score_scale = 3,
    icu_hd_ap = 3,
    pathogen_combined_types = 3,
    infection_types = 3
  ),
  # For pathogens[4] == "pse_car"
  list(
    age_new = 3,
    country_region = 3,
    country_income = 3,
    
    icu_hd_ap = 3,
    pathogen_combined_types = 3,
    infection_types = 3
  ),
  # For pathogens[5] == "entc_van"
  list(
    comorbidities_CCI = 3
  ),
  # For pathogens[6] == "sa_meth"
  list(
    severity_score_scale = 3
  ),
  # For pathogens[7] == "mdr_gnb"
  list(
    age_new = 3,
    country_region = 3,
    country_income = 3,
    comorbidities_CCI = 3,
    severity_score_scale = 3,
    icu_hd_ap = 3,
    pathogen_combined_types = 3,
    infection_types = 3
  ),
  # For pathogens[8] == "mdr"
  list(
    age_new = 3,
    country_region = 3,
    country_income = 3,
    comorbidities_CCI = 3,
    severity_score_scale = 3,
    icu_hd_ap = 3,
    pathogen_combined_types = 3,
    infection_types = 3
  ),
  # For pathogens[9] == "amr"
  list(
    age_new = 3,
    country_region = 3,
    country_income = 3,
    comorbidities_CCI = 3,
    severity_score_scale = 3,
    icu_hd_ap = 3,
    pathogen_combined_types = 3,
    infection_types = 3
  )
)


# Initialize a data.frame to store results
results_df <- data.frame(
  pathogen = character(),
  Hazard_Ratio = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each pathogen subset
for (i in seq_along(pathogens)) {
  var_name <- pathogens[i]
  df_used <- df_list[[i]]
  tvc_specs <- tvc_list[[i]]
  
  #Special handling for VRE
  if (var_name == "entc_van") {
    df_used$infection_types <- droplevels(factor(df_used$infection_types))
    df_used$infection_types <- relevel(df_used$infection_types, ref = "Hospital-acquired BSI")
  }
  
  # Build the model formula
  formula_text <- paste0(
    "Surv(time, event) ~ age_new + sex + country_region + country_income + hpd_admreason + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + infection_types + ", var_name
  )
  
  # Fit the Royston–Parmar model
  model <- stpm2(
    formula = as.formula(formula_text),
    data = df_used,
    df = 5,
    tvc = tvc_specs
  )
  
  # Find the “resistant” coefficient
  coef_names <- names(coef(model))
  idx <- grep(paste0("^", var_name, ".*resistant$"), coef_names, ignore.case = TRUE)
  
  beta <- coef(model)[idx]
  se <- sqrt(vcov(model)[idx, idx])
  
  # HR, CI and raw p-value
  HR <- exp(beta)
  CI <- exp(beta + c(-1, 1) * 1.96 * se)
  z_stat <- beta / se
  p_raw <- 2 * (1 - pnorm(abs(z_stat)))
  
  # Format p-value
  p_fmt <- ifelse(
    p_raw < 0.001, "<0.001",
    ifelse(p_raw < 0.01, sprintf("%.3f", p_raw),
           sprintf("%.2f", p_raw))
  )
  
  
  # Append to results
  results_df <- bind_rows(results_df, 
                          tibble::tibble(
                            pathogen = var_name,
                            counts = nrow(df_used),
                            Hazard_Ratio = HR,
                            Lower_CI = CI[1],
                            Upper_CI = CI[2]
                          )
  )
}

# ---------------
# Incidence among non-exposed
calculate_incidence_rate <- function(df, var) {
  # find all levels ending in "-susceptible"
  lvl_susc <- grep("-susceptible$", unique(df[[var]]), value = TRUE)
  if (length(lvl_susc) == 0) return(NA_real_)
  
  df0 <- df[df[[var]] %in% lvl_susc, ]
  # proportion of events in the susceptible group
  mean(df0$event == 1)
}

I0_df <- lapply(pathogens, function(var) {
  data.frame(
    pathogen = var,
    I0_rate  = calculate_incidence_rate(df_list[[which(pathogens == var)]], var),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Prevalence in the entire population 
pre_df <- lapply(pathogens, function(var) {
  # count levels ending in "-resistant"
  n_res <- sum(grepl("-resistant$", pop[[var]]),   na.rm = TRUE)
  data.frame(
    pathogen = var,
    prevalence = n_res / nrow(pop),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Combine
final_df <- results_df %>%
  left_join(I0_df, by = "pathogen") %>%
  left_join(pre_df, by = "pathogen")

# -----------------
# pathogen, Hazard_Ratio, Lower_CI, Upper_CI, I0_rate, prevalence
att <- final_df %>%
  mutate(
    HR_CI = sprintf("%.2f [%.2f, %.2f]", Hazard_Ratio, Lower_CI, Upper_CI),
    Incidence = sprintf("%.2f%%", I0_rate     * 100),
    Prevalence = sprintf("%.2f%%", prevalence * 100),
    # AR, AF, PAF and their CIs on the raw scale
    AR = I0_rate * (Hazard_Ratio - 1),
    AR_lo = I0_rate * (Lower_CI - 1),
    AR_hi = I0_rate * (Upper_CI - 1),
    AF = (Hazard_Ratio - 1) / Hazard_Ratio,
    AF_lo = (Lower_CI - 1) / Lower_CI,
    AF_hi = (Upper_CI - 1) / Upper_CI,
    PAF = (prevalence * (Hazard_Ratio - 1)) /
      (1 + prevalence * (Hazard_Ratio - 1)),
    PAF_lo = (prevalence * (Lower_CI - 1)) /
      (1 + prevalence * (Lower_CI - 1)),
    PAF_hi = (prevalence * (Upper_CI - 1)) /
      (1 + prevalence * (Upper_CI - 1)),

    AR_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", AR  * 100, AR_lo * 100, AR_hi * 100),
    AF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", AF  * 100, AF_lo * 100, AF_hi * 100),
    PAF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", PAF * 100, PAF_lo * 100, PAF_hi * 100)
  ) 

# For table
att$Group <- "Overall"
att$num <- paste0("N = ", att$counts)

# Save
saveRDS(att, "data/att28_paf_overall_conditional.RData")
###
