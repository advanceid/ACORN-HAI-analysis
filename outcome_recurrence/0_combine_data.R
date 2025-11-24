# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr,
                 dplyr,
                 tidyr,
                 gt,
                 gtsummary,
                 htmltools, 
                 openxlsx)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
baseline <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
episode <- readRDS("data/clean_data_RData/ast_all.RData")

# Filter episodes with both values 1 and 2
filter_multi_episodes <- function(data, col_name) {
  valid_ids <- data %>%
    group_by(recordid) %>%
    filter(any(.data[[col_name]] == 1) & any(.data[[col_name]] == 2)) %>%
    pull(recordid) %>%
    unique()
  
  data %>%
    filter(recordid %in% valid_ids & .data[[col_name]] %in% c(1,2))
}

id_vap <- baseline$recordid[grepl("VAP", baseline$infection_types)]
id_bsi <- baseline$recordid[grepl("BSI$", baseline$infection_types)]

# Filter episodes for VAP and BSI
vap_episode <- episode %>% filter(recordid %in% id_vap)
bsi_episode <- episode %>% filter(recordid %in% id_bsi)

# Episodes with both 1 and 2 for VAP and BSI
df_vap_recur <- filter_multi_episodes(vap_episode, "f07m_vapc")
df_bsi_recur <- filter_multi_episodes(bsi_episode, "f07m_bsic")

# Only 1
get_remain_episode <- function(episode_data, recur_data) {
  episode_data %>%
    filter(!recordid %in% recur_data$recordid) %>%
    group_by(recordid) %>%
    slice_min(order_by = as.Date(spec_date), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(recordid, spec_date)
}

df_vap_remain <- get_remain_episode(vap_episode, df_vap_recur)
df_bsi_remain <- get_remain_episode(bsi_episode, df_bsi_recur)

# ------------------
# BSI persistence
matched_org <- df_bsi_recur %>%
  filter(f07m_bsic %in% c(1,2)) %>%
  group_by(recordid, org_names_all) %>%
  summarise(has_1 = any(f07m_bsic == 1),
            has_2 = any(f07m_bsic == 2),
            .groups = "drop") %>%
  filter(has_1 & has_2)

# ris_
antibiotic_cols <- grep("^ris_", colnames(df_bsi_recur), value = TRUE)

# recordid, org_names_all, ris_
df_bsi_check <- df_bsi_recur %>%
  semi_join(matched_org, by = c("recordid", "org_names_all")) %>%
  filter(f07m_bsic %in% c(1,2)) %>%
  group_by(recordid, org_names_all) %>%
  summarise(across(all_of(antibiotic_cols), ~ n_distinct(.x) == 1),
            .groups = "drop")

consistent_antibiotic <- df_bsi_check %>%
  rowwise() %>%
  filter(all(c_across(all_of(antibiotic_cols)))) %>%
  ungroup() %>%
  select(recordid, org_names_all)

df_bsi_persist <- df_bsi_recur %>%
  semi_join(consistent_antibiotic, by = c("recordid", "org_names_all")) %>%
  filter(f07m_bsic %in% c(1,2))
# ------------
# Function to calculate duration between spec_dates (generalized)
calculate_duration <- function(data, id_col, event_col, date_col, prefix) {
  data %>%
    group_by(across(all_of(c(id_col, event_col)))) %>%
    summarise(first_date = min(as.Date(.data[[date_col]])), .groups = "drop") %>%
    pivot_wider(names_from = all_of(event_col), 
                values_from = first_date, 
                names_prefix = "spec_date_") %>%
    mutate(!!paste0("duration_", prefix) := as.numeric(.data[["spec_date_2"]] - 
                                                         .data[["spec_date_1"]])) %>%
    select(all_of(id_col), starts_with("spec_date_"), starts_with("duration_"))
}

# 
bsi_persist_duration <- calculate_duration(df_bsi_persist, "recordid", "f07m_bsic", "spec_date", "bsi_persist")
vap_recur_duration <- calculate_duration(df_vap_recur, "recordid", "f07m_vapc", "spec_date", "vap_recur")
bsi_recur_duration <- calculate_duration(df_bsi_recur, "recordid", "f07m_bsic", "spec_date", "bsi_recur")

# Add baseline
process_recurrence <- function(baseline_df, ids, duration_df, indicator_name) {
  baseline_df %>%
    filter(recordid %in% ids) %>%
    mutate("{{indicator_name}}" := ifelse(recordid %in% duration_df$recordid, 1, 0)) %>%
    left_join(duration_df, by = "recordid")
}

# 
vap_recur <- process_recurrence(baseline, id_vap, vap_recur_duration, recurrence)
bsi_recur <- process_recurrence(baseline, id_bsi, bsi_recur_duration, recurrence)


# Function to fill missing spec_date_1 with spec_date from remain episodes (only 1)
fill_spec_date <- function(recur_df, remain_df, spec_date_col = "spec_date") {
  recur_df %>%
    left_join(remain_df %>% select(recordid, remain_spec_date = !!sym(spec_date_col)), by = "recordid") %>%
    mutate(spec_date_1 = if_else(is.na(spec_date_1), remain_spec_date, spec_date_1)) %>%
    select(-remain_spec_date)
}

# Delete spec_date NA
vap_recur_final <- fill_spec_date(vap_recur, df_vap_remain) %>%
  filter(!is.na(spec_date_1))

bsi_recur_final <- fill_spec_date(bsi_recur, df_bsi_remain) %>%
  filter(!is.na(spec_date_1))

###
# 1 (recurrence), 2 (dead)
assign_event <- function(df) {
  df %>%
    mutate(event = case_when(
      !is.na(mortality_date) & mortality_date >= spec_date_1 & is.na(spec_date_2) ~ 2,
      !is.na(mortality_date) & !is.na(spec_date_2) & mortality_date >= spec_date_2 ~ 1,
      is.na(mortality_date) & !is.na(spec_date_2) ~ 1,
      is.na(mortality_date) & is.na(spec_date_2) ~ 0,
      TRUE ~ NA_real_
    ))
}

vap_recur_final <- assign_event(vap_recur_final)
bsi_recur_final <- assign_event(bsi_recur_final)

# Time 
process_time_recur <- function(df, duration_col) {
  df %>%
    mutate(
      dead_spec = as.numeric(mortality_date - spec_date_1),
      time_ori = case_when(
        event == 1 ~ .data[[duration_col]],
        event == 2 ~ dead_spec,
        event == 0 ~ 90
      ),
      event = if_else(time_ori > 90, 0, event),
      time = if_else(event == 0, 90, time_ori)
    ) %>% 
    filter(!is.na(event))
}

vap_recur_final <- process_time_recur(vap_recur_final, "duration_vap_recur")
bsi_recur_final <- process_time_recur(bsi_recur_final, "duration_bsi_recur")

# BSI persistence
process_bsi_persistence <- function(df, persist_ids) {
  df %>%
    mutate(
      time_new = if_else(
        event == 1 & !recordid %in% persist_ids,
        pmin(dead_spec, 90),
        time
      ),
      event_new = if_else(
        event == 1 & (!recordid %in% persist_ids) & dead_spec <= 90,
        2,
        event
      ),
      event_new = replace_na(event_new, 0),
      time_new = if_else(event_new == 0, 90, time_new)
    )
}

#
bsi_persist_final <- process_bsi_persistence(
  bsi_recur_final, 
  bsi_persist_duration$recordid
)

#
all_data <- list(vap_recur_final, bsi_recur_final, bsi_persist_final)

# Save
saveRDS(all_data, "data/rp_data.RData")
###
