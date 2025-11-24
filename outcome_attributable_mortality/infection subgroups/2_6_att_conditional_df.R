# Clear
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    survival,
    rstpm2,
    gt,
    extrafont
  )
})

# Set working directory and load fonts
setwd("./")
loadfonts()

# Load data
df_list <- readRDS("data/att_first28_infection.RData")

# Candidate spline degrees of freedom
df_sp <- 1:3

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

# Prepare a list to collect all AIC data.frames
all_aic <- list()

# Loop over each pathogen
for (var_name in pathogens) {
  
  # Find the matching dataset
  df_ori <- df_list[[ which(pathogens == var_name) ]]
  
  # Loop over each income stratum
  for (inc_pref in inf_levels_prefixed) {
    # Strip numeric prefix
    inc <- sub("^\\d+\\s+", "", inc_pref)
    
    # Skip if stratum not present
    if (!inc %in% df_ori$infection_types) next
    
    # Subset to this income level
    df_used <- df_ori %>%
      filter(infection_types == inc)
    
    # Count total obs and resistant cases
    n_obs    <- nrow(df_used)
    n_resist <- sum(
      grepl("resistant$", df_used[[var_name]]),
      na.rm = TRUE
    )
    
    # Skip under-powered strata
    if (n_obs < 100)  next
    if (n_resist < 30)   next
    
    # Compute AIC for each candidate df
    aic_vals <- vapply(df_sp, function(dfi) {
      # Build formula text
      ftxt <- paste0(
        "Surv(time, event) ~ age_new + sex + country_region + country_income + hpd_admreason + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + ", var_name
      )
      mod <- stpm2(
        formula = as.formula(ftxt),
        data = df_used,
        df = dfi
      )
      AIC(mod)
    }, numeric(1))
    
    # Assemble and sort
    aic_df <- tibble::tibble(
      DF = df_sp,
      AIC = aic_vals,
      pathogen = var_name,
      income_level = inc
    ) %>%
      arrange(AIC)
    
    # Store 
    all_aic[[ paste(var_name, inc, sep = " | ") ]] <- aic_df
  }
}

# Combine all 
aic_master <- bind_rows(all_aic, .id = "stratum")

