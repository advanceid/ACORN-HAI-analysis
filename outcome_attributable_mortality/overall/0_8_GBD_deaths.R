# Clear 
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(openxlsx, 
                 dplyr, 
                 magrittr,
                 tidyr,
                 stringr)
})

# Set working directory
wd <- "./"
setwd(wd)

# Load data without colnames
df_raw <- read.xlsx("data/IHME_GBD_2021_MORTALITY_1990_2021_SR_TABLE_1_Y2024M04D03.xlsx", 
                    sheet = 1, colNames = FALSE)

# Remove the first row
df_tmp <- df_raw[-1, ]

# Use the second row as colnames
colnames(df_tmp) <- as.character(unlist(df_tmp[1, ]))

# Remove the header row (now colnames)
df <- df_tmp[-1, ]

# Keep only "Country" rows
df <- df %>% filter(location_type == "Country")

# Define vector of Asian countries
asian_countries <- c(
  "Afghanistan","Armenia","Azerbaijan","Bahrain","Bangladesh","Bhutan","Brunei Darussalam",
  "Cambodia","China","Cyprus","Democratic People's Republic of Korea","Georgia","India",
  "Indonesia","Iran (Islamic Republic of)","Iraq","Israel","Japan","Jordan","Kazakhstan",
  "Kuwait","Kyrgyzstan","Lao People's Democratic Republic","Lebanon","Malaysia","Maldives",
  "Mongolia","Myanmar","Nepal","Oman","Pakistan","Palestine","Philippines","Qatar",
  "Republic of Korea","Saudi Arabia","Singapore","Sri Lanka","Syrian Arab Republic",
  "Tajikistan","Thailand","Timor-Leste","Türkiye","Turkmenistan","United Arab Emirates",
  "Uzbekistan","Viet Nam","Yemen"
)

# Filter to Asian countries
df_asia <- df %>% filter(location_name %in% asian_countries)

# Add income level (World Bank FY2025 classification)
df_asia <- df_asia %>%
  mutate(income_level = case_when(
    # --- Low income ---
    location_name %in% c("Afghanistan","Yemen",
                         "Democratic People's Republic of Korea") ~ "Low income",
    
    # --- Lower middle income ---
    location_name %in% c("Bangladesh","Bhutan","Cambodia","India","Kyrgyzstan",
                         "Lao People's Democratic Republic","Myanmar","Nepal",
                         "Pakistan","Philippines","Tajikistan","Uzbekistan","Viet Nam",
                         "Palestine") ~ "Lower middle income",
    
    # --- Upper middle income ---
    location_name %in% c("Armenia","Azerbaijan","China","Georgia","Indonesia",
                         "Iran (Islamic Republic of)","Iraq","Jordan","Kazakhstan",
                         "Lebanon","Malaysia","Maldives","Mongolia","Sri Lanka",
                         "Syrian Arab Republic","Thailand","Timor-Leste",
                         "Turkmenistan","Türkiye") ~ "Upper middle income",
    
    # --- High income ---
    location_name %in% c("Bahrain","Brunei Darussalam","Cyprus","Israel","Japan",
                         "Kuwait","Oman","Qatar","Republic of Korea","Saudi Arabia",
                         "Singapore","United Arab Emirates") ~ "High income",
    
    TRUE ~ "Unclassified"
  ))

# Filter VAP, BSI
df_used <- df_asia %>%
  filter(cause_name %in% c("Lower respiratory infections", 
                           "Maternal sepsis and other maternal infections", 
                           "Neonatal sepsis and other neonatal infections")) %>%
  select(location_name, cause_name, `2019 (All-Age Deaths)`, income_level)

# Parse "2019 (All-Age Deaths)" into estimate / ci_low / ci_high (numeric)
df_used <- df_used %>%
  mutate(deaths_str = as.character(`2019 (All-Age Deaths)`)) %>%
  extract(
    col = deaths_str,
    into = c("estimate", "ci_low", "ci_high"),
    regex = "^\\s*([0-9.,]+)\\s*\\(([-0-9.,]+)\\s*-\\s*([-0-9.,]+)\\)\\s*$",
    remove = FALSE
  ) %>%
  mutate(
    estimate = as.numeric(str_replace_all(estimate, ",", "")) / 100000,
    ci_low   = as.numeric(str_replace_all(ci_low,   ",", "")) / 100000,
    ci_high  = as.numeric(str_replace_all(ci_high,  ",", "")) / 100000
  ) %>%
  select(-`2019 (All-Age Deaths)`, -deaths_str)

# Map causes to syndromes (LRI vs BSI)
df_used <- df_used %>%
  mutate(syndrome = case_when(
    cause_name == "Lower respiratory infections" ~ "LRI",
    cause_name %in% c("Maternal sepsis and other maternal infections",
                      "Neonatal sepsis and other neonatal infections") ~ "BSI",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(syndrome))

# Aggregate by income level × syndrome
agg_income <- df_used %>%
  group_by(income_level, syndrome) %>%
  summarise(
    deaths_est  = sum(estimate, na.rm = TRUE),
    deaths_low  = sum(ci_low,   na.rm = TRUE),
    deaths_high = sum(ci_high,  na.rm = TRUE),
    .groups = "drop"
  )

# Asia overall × syndrome
agg_asia_syndrome <- df_used %>%
  group_by(syndrome) %>%
  summarise(
    deaths_est  = sum(estimate, na.rm = TRUE),
    deaths_low  = sum(ci_low,   na.rm = TRUE),
    deaths_high = sum(ci_high,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(income_level = "Asia overall") %>%
  select(income_level, everything())

# Asia overall × all infection-related (LRI + BSI)
agg_asia_total <- df_used %>%
  summarise(
    deaths_est  = sum(estimate, na.rm = TRUE),
    deaths_low  = sum(ci_low,   na.rm = TRUE),
    deaths_high = sum(ci_high,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(income_level = "Asia overall",
         syndrome = "All infections") %>%
  select(income_level, syndrome, everything())

# Combine everything
agg_all <- bind_rows(agg_income, agg_asia_syndrome, agg_asia_total) %>%
  arrange(factor(income_level,
                 levels = c("Asia overall","High income","Upper middle income",
                            "Lower middle income","Low income","Unclassified")),
          syndrome)

# Save
saveRDS(agg_all, "data/GBD_death.RDS")
#

