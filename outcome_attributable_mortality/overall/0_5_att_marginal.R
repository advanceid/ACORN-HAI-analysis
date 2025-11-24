# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 coxphw,
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
  df_used  <- df_list[[i]]
  
  # Survival formula 
  surv_formula <- as.formula(paste("Surv(time, event) ~", var_name))
  
  # Fit weighted Cox model (average hazard ratio)
  model <- coxphw(
    formula = surv_formula,
    data = df_used,
    template = "AHR",   
    caseweights = df_used$wt,  
    robust = TRUE
  )
  
  # Extract coefficients and confidence intervals
  s <- summary(model)
  hr <- exp(s$coefficients)[1]  
  lo <- s$ci.lower[1]  
  up <- s$ci.upper[1] 
  
  # Append to results
  results_df <- bind_rows(results_df, 
                          tibble::tibble(
                            pathogen = var_name,
                            counts = nrow(df_used),
                            Hazard_Ratio = hr,
                            Lower_CI = lo,
                            Upper_CI = up
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
    pathogen   = var,
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
    Incidence = sprintf("%.2f%%", I0_rate * 100),
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
    # format them as "% [low, high]"
    AR_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", AR * 100, AR_lo * 100, AR_hi * 100),
    AF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", AF * 100, AF_lo * 100, AF_hi * 100),
    PAF_CI = sprintf("%.2f%% [%.2f%%, %.2f%%]", PAF * 100, PAF_lo * 100, PAF_hi * 100)
  ) 

# For table
att$Group <- "Overall"
att$num <- paste0("N = ", att$counts)

# Save
saveRDS(att, "data/att28_paf_overall_marginal.RData")
###


