# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr,
    dplyr,
    purrr
  )
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")

# Dates
df$inf_onset <- as.Date(df$inf_onset)
df$ho_discharge_date <- as.Date(df$ho_discharge_date)
df$mortality_date <- as.Date(df$mortality_date)

# Build time & event relative to infection onset
df_used_ori <- df %>%
  # remove any rows with missing infection onset
  dplyr::filter(!is.na(inf_onset)) %>%
  # ensure transfer or death not before onset
  dplyr::filter(
    is.na(ho_discharge_date) | ho_discharge_date >= inf_onset,
    is.na(mortality_date) | mortality_date >= inf_onset
  ) %>%
  # raw intervals
  dplyr::mutate(
    t_out = as.numeric(ho_discharge_date - inf_onset),
    t_death = as.numeric(mortality_date - inf_onset)
  ) %>%
  # convert any 0 days to 1 day
  dplyr::mutate(
    t_out   = dplyr::if_else(!is.na(t_out) & t_out   == 0, 1, t_out),
    t_death = dplyr::if_else(!is.na(t_death) & t_death == 0, 1, t_death)
  ) %>%
  # derive min time and event code: 1 = discharge, 2 = death, 0 = censored
  dplyr::mutate(
    time = pmin(t_out, t_death, na.rm = TRUE),
    event = dplyr::case_when(
      !is.na(t_out) & (is.na(t_death) | t_out < t_death) ~ 1L,
      !is.na(t_death) & (is.na(t_out) | t_death <= t_out) ~ 2L,
      TRUE ~ 0L
    )
  )

# Cut-off helper
process_cutoff <- function(data, t) {
  data %>%
    dplyr::mutate(
      time_new  = dplyr::if_else(is.na(time) | time > t, t,  time),
      event_new = dplyr::if_else(is.na(time) | time > t, 0L, event)
    )
}

# Cut-offs
cutoffs <- c(14, 21, 28, 60, 90)

# List of datasets for each cutoff
df_used <- purrr::map(cutoffs, ~ process_cutoff(df_used_ori, .x))
names(df_used) <- paste0("t", cutoffs)

# ---------------------------
# Subset helper
filter_non_na <- function(df, column_name,
                          trim = TRUE,
                          trim_p = c(0.01, 0.99)) {
  # Check column exists
  if (!column_name %in% colnames(df)) {
    stop(paste("Column", column_name, "does not exist in the data frame."))
  }
  
  # Filter out NA in target column and create ast_ris
  df <- df %>%
    dplyr::filter(!is.na(.data[[column_name]])) %>%
    dplyr::mutate(
      ast_ris = ifelse(grepl("-resistant$", .data[[column_name]]), 1, 0),
      ast_ris = factor(ast_ris, levels = c(0, 1))
    )
  
  # Relevel icu_hd_ap to "No" if present
  if ("icu_hd_ap" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(icu_hd_ap = stats::relevel(factor(icu_hd_ap), ref = "No"))
  }
  
  # Drop unused levels for pathogen_combined_types if present
  if ("pathogen_combined_types" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(pathogen_combined_types = droplevels(factor(pathogen_combined_types)))
  }
  
  df
}

# Pathogen columns and desired names
pathogens    <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth")
subset_names <- c("car_aci", "thir_ceph_ent", "car_ent", "car_pse", "van_entc", "sa_meth")

# Build list-of-lists of subsets for each cutoff
all_subsets <- purrr::map(df_used, function(dat) {
  purrr::map(pathogens, ~ filter_non_na(dat, .x)) %>%
    purrr::set_names(subset_names)
})

# ---------------------------
# Post-process: ONLY for entc_van subset (named "van_entc"),
# drop unused levels of infection_types
all_subsets <- purrr::map(all_subsets, function(lst) {
  if (!is.null(lst$van_entc) && "infection_types" %in% names(lst$van_entc)) {
    lst$van_entc <- lst$van_entc %>%
      dplyr::mutate(infection_types = factor(infection_types)) %>%
      dplyr::mutate(infection_types = droplevels(infection_types))
  }
  lst
})

# Save
saveRDS(all_subsets, file = "data/los_lists.RData")



# Discharge
los_summary <- df_used_ori %>%
  filter(!is.na(t_out)) %>%
  summarise(
    median_LOS = median(t_out, na.rm = TRUE),
    IQR_LOS    = IQR(t_out, na.rm = TRUE),
    Q1_LOS     = quantile(t_out, 0.25, na.rm = TRUE),
    Q3_LOS     = quantile(t_out, 0.75, na.rm = TRUE),
    n          = n()
  )

los_summary

