# Clear 
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    magrittr,
    survival,
    purrr,
    rstpm2,
    tibble,
    scales
  )
})

# Set working directory
setwd("./")

# Load data
pop <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
df_list <- readRDS("data/att_first28_infection.RData")

# Define pathogen variables
pathogens <- c(
  "aci_car",
  "ent_thir",
  "ent_car",
  "pse_car",
  "entc_van",
  "sa_meth",
  "mdr_gnb"
)

# Define infection strata 
inf_levels_prefixed <- c(
  "1 VAP",
  "2 Hospital-acquired BSI",
  "3 Healthcare-associated BSI"
)

# Define nested tvc settings: first by pathogen, then by infection level
tvc_list <- list(
  aci_car = list(
    "VAP" = list(age_new = 2, sex = 2, pathogen_combined_types = 2),
    "Hospital-acquired BSI" = list(comorbidities_CCI = 2, icu_hd_ap = 2),
    "Healthcare-associated BSI" = NULL
  ),
  ent_thir = list(
    "VAP" = list(country_region = 2),
    "Hospital-acquired BSI" = list(severity_score_scale = 2),
    "Healthcare-associated BSI" = list(age_new = 2, country_income = 2)
  ),
  ent_car = list(
    "VAP" = list(country_region = 2),
    "Hospital-acquired BSI" = list(age_new = 2, country_region = 2, country_income = 2, severity_score_scale = 2),
    "Healthcare-associated BSI" = list(age_new = 2, sex = 2, country_income = 2
    )
  ),
  pse_car = list(
    "VAP" = list(country_region = 2),
    "Hospital-acquired BSI" = list(age_new = 2, sex = 2, country_region = 2, country_income = 2, severity_score_scale = 2),
    "Healthcare-associated BSI" = NULL
  ),
  entc_van = list(
    "VAP" = NULL,
    "Hospital-acquired BSI" = list(comorbidities_CCI = 2, severity_score_scale = 2),
    "Healthcare-associated BSI" = NULL
  ),
  sa_meth = list(
    "VAP" = NULL,
    "Hospital-acquired BSI" = list(severity_score_scale = 2),
    "Healthcare-associated BSI" = list(age_new = 2, country_region = 2)
  ),
  mdr_gnb = list(
    "VAP" = list(age_new = 2, country_region = 2),
    "Hospital-acquired BSI" = list(age_new = 2, country_region = 2, country_income = 2, 
                                   comorbidities_CCI = 2, severity_score_scale = 2, icu_hd_ap = 2),
    "Healthcare-associated BSI" = list(age_new = 2, country_income = 2, severity_score_scale = 2)
  )
)

# Initialize results tibble
results_df <- tibble(
  pathogen = character(),
  infection = character(),
  counts = integer(),
  Hazard_Ratio = double(),
  Lower_CI = double(),
  Upper_CI = double()
)

# Loop over each pathogen and each income stratum
for (var_name in pathogens) {
  # retrieve the matching dataset
  df_ori <- df_list[[ which(pathogens == var_name) ]]
  
  for (inc_pref in inf_levels_prefixed) {
    # strip numeric prefix
    inc <- sub("^\\d+\\s+", "", inc_pref)
    # skip if this infection stratum is not present
    if (!inc %in% df_ori$infection_types) next
    
    # subset to this infection level
    df_used <- df_ori %>% filter(infection_types == inc)
    
    # count observations and resistant cases
    n_obs    <- nrow(df_used)
    n_resist <- sum(grepl("resistant$", df_used[[var_name]]), na.rm = TRUE)
    
    # skip under‐powered strata
    if (n_obs < 100) next
    if (n_resist < 30) next
    
    # pick the tvc spec for this subgroup
    tvc_specs <- tvc_list[[var_name]][[inc]]
    
    # build the model formula
    formula_text <- paste0(
      "Surv(time, event) ~ age_new + sex + country_region + country_income + hpd_admreason + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + ", var_name
    )
    
    # fit the Royston–Parmar model
    model <- stpm2(
      formula = as.formula(formula_text),
      data = df_used,
      df = 2,
      tvc = tvc_specs
    )
    
    # extract the coefficient and SE for the "resistant" level
    coef_names <- names(coef(model))
    idx        <- grep(
      paste0("^", var_name, ".*resistant$"),
      coef_names,
      ignore.case = TRUE
    )
    beta <- coef(model)[idx]
    se   <- sqrt(vcov(model)[idx, idx])
    
    # compute HR and 95% CI
    HR <- exp(beta)
    CI <- exp(beta + c(-1, 1) * 1.96 * se)
    
    # append to results
    results_df <- bind_rows(
      results_df,
      tibble(
        pathogen = var_name,
        infection = inc,
        counts = n_obs,
        Hazard_Ratio = HR,
        Lower_CI = CI[1],
        Upper_CI = CI[2]
      )
    )
  }
}

# ---------------
# Incidence among non-exposed
I0_df <- map_dfr(pathogens, function(var) {
  map_dfr(inf_levels_prefixed, function(inc_pref) {
    # strip prefix
    inc <- sub("^\\d+\\s+", "", inc_pref)
    # subset the df for this pathogen & infection
    df_var <- df_list[[which(pathogens == var)]] %>%
      filter(infection_types == inc)
    # compute susceptible‐group event rate
    lvl_susc <- grep("-susceptible$", unique(df_var[[var]]), value = TRUE)
    I0_rate <- if (length(lvl_susc) == 0) {
      NA_real_
    } else {
      df_susc <- df_var %>% filter(.data[[var]] %in% lvl_susc)
      mean(df_susc$event == 1, na.rm = TRUE)
    }
    # return one tibble row
    tibble(
      pathogen = var,
      infection_pref = inc_pref,
      I0_rate = I0_rate
    )
  })
}) %>%
  mutate(infection = sub("^\\d+\\s+", "", infection_pref)) %>%
  select(-infection_pref)


# Prevalence of “-resistant”, by pathogen & income
pre_df <- map_dfr(pathogens, function(var) {
  map_dfr(inf_levels_prefixed, function(inc_pref) {
    
    inc <- sub("^\\d+\\s+", "", inc_pref)
    
    pop_sub <- pop %>% filter(infection_types == inc)
    n_res <- sum(grepl("-resistant$", pop_sub[[var]]), na.rm = TRUE)
    
    prevalence <- n_res / nrow(pop_sub)
    
    tibble(
      pathogen = var,
      infection_pref = inc_pref,
      prevalence = prevalence
    )
  })
}) %>%
  mutate(infection = sub("^\\d+\\s+", "", infection_pref)) %>%
  select(-infection_pref)

# Combine
final_df <- results_df %>%    
  left_join(I0_df,  by = c("pathogen", "infection")) %>%
  left_join(pre_df, by = c("pathogen", "infection"))

# ------------
# pathogen, Hazard_Ratio, Lower_CI, Upper_CI, I0_rate, prevalence
att <- final_df %>%
  mutate(
    HR_CI = sprintf("%.2f [%.2f, %.2f]", Hazard_Ratio, Lower_CI, Upper_CI),
    Incidence = sprintf("%.2f%%", I0_rate * 100),
    Prevalence = sprintf("%.2f%%", prevalence * 100),
    # AR, AF, PAF and CIs on the raw scale
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
    # Combine
    AR_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", AR * 100, AR_lo * 100, AR_hi * 100),
    AF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", AF * 100, AF_lo * 100, AF_hi * 100),
    PAF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", PAF * 100, PAF_lo * 100, PAF_hi * 100)
  ) 

# For table
att$Group <- "Infection syndromes"
att$num <- paste0(att$infection, " (", att$counts, ")")

# Save
saveRDS(att, "data/att28_paf_infection_conditional.RData")
###