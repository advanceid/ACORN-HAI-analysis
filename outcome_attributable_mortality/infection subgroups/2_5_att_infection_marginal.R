# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, 
                 magrittr,
                 purrr,
                 coxphw, 
                 scales, 
                 tibble,
                 openxlsx)
})

# Set working directory
wd <- "./"
setwd(wd)

# Read data
pop <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
df_list <- readRDS("data/att_first28_infection.RData")

# Define pathogen variables
pathogens <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth", "mdr_gnb")

# Define fixed infection levels
infection_levels_prefixed <- c(
  "1 VAP",
  "2 Hospital-acquired BSI",
  "3 Healthcare-associated BSI"
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

# Loop over each pathogen and each prefixed infection level
for (var_name in pathogens) {
  df_ori <- df_list[[which(pathogens == var_name)]]
  
  for (inc_pref in infection_levels_prefixed) {
    # strip numeric prefix to match the column values
    inc <- sub("^\\d+\\s+", "", inc_pref)
    # Skip if this income stratum is absent
    if (!inc %in% df_ori$infection_types) next
    
    df_used <- df_ori %>% filter(infection_types == inc)
    
    # Count total observations and “-resistant” cases
    n_obs    <- nrow(df_used)
    n_resist <- sum(
      grepl("resistant$", df_used[[var_name]]),
      na.rm = TRUE
    )
    
    # Skip small or under‐powered strata
    if (n_obs < 100)   next
    if (n_resist < 30)   next
    
    # Build survival formula
    surv_formula <- as.formula(paste("Surv(time, event) ~", var_name))
    
    # Fit weighted Cox model
    model <- coxphw(
      formula = surv_formula,
      data = df_used,
      template = "AHR",
      caseweights = df_used$wt,
      robust = TRUE
    )
    
    # Extract HR and 95% CI
    s  <- summary(model)
    hr <- exp(s$coefficients)[1]
    lo <- s$ci.lower[1]
    up <- s$ci.upper[1]
    
    # Append to results
    results_df <- bind_rows(
      results_df,
      tibble(
        pathogen = var_name,
        infection = inc_pref,
        counts = n_obs,
        Hazard_Ratio = hr,
        Lower_CI = lo,
        Upper_CI = up
      )
    )
  }
}

# Enforce factor ordering for both pathogen and infection, then sort
results_df <- results_df %>%
  mutate(
    pathogen = factor(pathogen, levels = pathogens),
    infection = factor(infection, levels = infection_levels_prefixed)
  ) %>%
  arrange(pathogen, infection) %>%
  # finally drop the numeric prefix for display
  mutate(infection = sub("^\\d+\\s+", "", as.character(infection)))

# ---------------
# Incidence among non-exposed
I0_df <- map_dfr(pathogens, function(var) {
  map_dfr(infection_levels_prefixed, function(inc_pref) {
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


# Prevalence of “-resistant”, by pathogen & infection
pre_df <- map_dfr(pathogens, function(var) {
  map_dfr(infection_levels_prefixed, function(inc_pref) {

    inc <- sub("^\\d+\\s+", "", inc_pref)

    pop_sub <- pop %>% filter(infection_types == inc)
    n_res <- sum(grepl("-resistant$", pop_sub[[var]]), na.rm = TRUE)

    prevalence <- n_res / nrow(pop_sub)

    tibble(
      pathogen = var,
      infection_pref = inc_pref,
      prevalence  = prevalence
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
    AR_CI  = sprintf("%.2f%% [%.2f%%, %.2f%%]", AR  * 100, AR_lo * 100, AR_hi * 100),
    AF_CI  = sprintf("%.2f%% [%.2f%%, %.2f%%]", AF  * 100, AF_lo * 100, AF_hi * 100),
    PAF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", PAF * 100, PAF_lo * 100, PAF_hi * 100)
  ) 

# For table
att$Group <- "Infection syndromes"
att$num <- paste0(att$infection, " (",att$counts, ")")

# Save
saveRDS(att, "data/att28_paf_infection_marginal.RData")
###
att_mHR <- att[,c(1,2,9)]

write.xlsx(att_mHR, "att28_mHR.xlsx", rowNames = FALSE)
