# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    openxlsx, dplyr, magrittr, purrr, lubridate, stringr,
    tidyverse, lme4, meta, ggplot2, tibble, grid, Cairo, tidyr
  )
})

# Paths & data
wd <- "./"
setwd(wd)

# Load data
df <- read.csv("data/ACORNHAISiteSummary-Wards_DATA_LABELS_2025-07-22_0938.csv")
study_period <- read.xlsx("data/Study duration.xlsx")
bed_estimate <- read.xlsx("data/ACORN-HAI beds.xlsx")
episode <- readRDS("data/clean_data_RData/episode_for_incidence.RData")

# Site mapping
df <- df %>% mutate(Site = substr(Site.Name, 1, 5))

all_siteids <- c(
  unique(as.character(episode$siteid)),
  unique(as.character(df$Site)),
  unique(as.character(study_period$Site))
)
all_siteids <- sort(unique(all_siteids[!is.na(all_siteids) & all_siteids != ""]))

site_map <- tibble(siteid = all_siteids) %>%
  arrange(siteid) %>%
  mutate(site_label = sprintf("Site %02d", row_number()))

df <- df %>% mutate(siteid = Site) %>% left_join(site_map, by = "siteid")
study_period <- study_period %>% mutate(siteid = Site) %>% left_join(site_map, by = "siteid")
episode <- episode %>% left_join(site_map, by = "siteid")

# ICU flag from episode
episode <- episode %>%
  mutate(icu_any = ifelse(!is.na(adm_ward_types_new) & adm_ward_types_new == "ICU", TRUE, FALSE))

# Keep only VAP + HAI-BSI
episode <- episode %>% filter(infection_types %in% c("VAP", "Hospital-acquired BSI"))


# Study period cleaning & gap removal
study_period <- study_period %>%
  drop_na() %>%
  mutate(
    Start.screening = as.Date(Start.screening, origin = "1899-12-30"),
    Recruitment.end = as.numeric(Recruitment.end),
    Recruitment.end = as.Date(Recruitment.end, origin = "1899-12-30"),
    Recruitment.end = ceiling_date(Recruitment.end, "month") - days(1),
    period = as.numeric(Recruitment.end - Start.screening + 1)
  )

episode_merged <- episode %>%
  left_join(study_period, by = c("siteid" = "Site")) %>%
  mutate(date_enrolment = as.Date(date_enrolment)) %>%
  arrange(siteid, date_enrolment) %>%
  group_by(siteid) %>%
  filter(hpd_adm_date >= Start.screening)

gap_info <- episode_merged %>%
  summarise(
    date_list = list(sort(date_enrolment)),
    Start.screening = first(Start.screening),
    Recruitment.end = first(Recruitment.end),
    .groups = "drop"
  ) %>%
  mutate(
    gap_periods = purrr::map(date_list, function(dates) {
      if (length(dates) < 2) return(tibble(start = as.Date(character()), end = as.Date(character())))
      gaps <- as.numeric(diff(dates))
      valid_idx <- which(gaps > 14)
      if (length(valid_idx) == 0) return(tibble(start = as.Date(character()), end = as.Date(character())))
      starts <- dates[-length(dates)][valid_idx] + days(1)
      ends   <- dates[-1][valid_idx] - days(1)
      tibble(start = starts, end = ends)
    })
  )

adjusted_periods <- gap_info %>%
  mutate(
    valid_periods = pmap(list(Start.screening, Recruitment.end, gap_periods),
                         function(start_screen, end_recruit, gaps) {
                           all_periods <- tibble(start = start_screen, end = end_recruit)
                           if (nrow(gaps) == 0) return(all_periods)
                           
                           for (i in seq_len(nrow(gaps))) {
                             gap <- gaps[i, ]
                             new_periods <- list()
                             for (j in seq_len(nrow(all_periods))) {
                               p <- all_periods[j, ]
                               if (gap$start > p$start && gap$end < p$end) {
                                 new_periods <- append(new_periods, list(
                                   tibble(start = p$start, end = gap$start - days(1)),
                                   tibble(start = gap$end + days(1), end = p$end)
                                 ))
                               } else if (gap$start <= p$start && gap$end < p$end && gap$end >= p$start) {
                                 new_periods <- append(new_periods, list(tibble(start = gap$end + days(1), end = p$end)))
                               } else if (gap$start > p$start && gap$start <= p$end && gap$end >= p$end) {
                                 new_periods <- append(new_periods, list(tibble(start = p$start, end = gap$start - days(1))))
                               } else if (gap$start <= p$start && gap$end >= p$end) {
                                 # remove
                               } else {
                                 new_periods <- append(new_periods, list(p))
                               }
                             }
                             all_periods <- bind_rows(new_periods)
                           }
                           all_periods
                         })
  )

final_valid_periods <- adjusted_periods %>%
  select(siteid, valid_periods) %>%
  unnest(valid_periods) %>%
  mutate(days = as.numeric(end - start + 1)) %>%
  filter(days >= 15)

valid_days_by_site <- final_valid_periods %>%
  group_by(siteid) %>%
  summarise(valid_period = sum(days), .groups = "drop")

study_period_updated <- study_period %>%
  rename(original_period = period) %>%
  left_join(valid_days_by_site, by = c("Site" = "siteid")) %>%
  mutate(
    period_valid_final = coalesce(valid_period, 0),
    reduced_days = original_period - period_valid_final
  )


# Ward categorization (for denominators)
df <- df %>%
  mutate(ward_cal = case_when(
    str_detect(Ward, regex("ICU|CCU|VCU|Care Unit", ignore_case = TRUE)) ~ "ICU",
    str_detect(Ward, regex("HDU|IMCU", ignore_case = TRUE)) ~ "IMCU",
    TRUE ~ "Other"
  ))


# Episode splits by syndrome & ICU
vap_episode      <- episode %>% filter(redcap_event_name == "vap_episode_arm_1")
bsi_episode      <- episode %>% filter(redcap_event_name == "bsi_episode_arm_1")

vap_episode_icu  <- vap_episode %>% filter(icu_any)
bsi_episode_all  <- bsi_episode
bsi_episode_icu  <- bsi_episode %>% filter(icu_any)

# Pathogen variable map
patho_vars <- c("aci_car","ent_thir","ent_car","pse_car","entc_van","sa_meth")
patho_map  <- c(aci_car="CRA", ent_thir="3GCRE", ent_car="CRE",
                pse_car="CRP", entc_van="VRE",  sa_meth="MRSA")
pathogen_list <- c("CRA","3GCRE","CRE","CRP","VRE","MRSA")

summarise_pathogen <- function(ep_df) {
  out <- lapply(patho_vars, function(v) {
    ep_df %>%
      filter(.data[[v]] == 1) %>%
      group_by(siteid, country_income) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(pathogen = patho_map[[v]])
  }) %>% bind_rows()
  
  if (nrow(out) == 0)
    return(tibble(siteid=character(), country_income=character()))
  pivot_wider(out, names_from = pathogen, values_from = n, values_fill = 0)
}

# Episode totals + pathogens
vap_summary <- vap_episode_icu %>%
  group_by(siteid, country_income) %>%
  summarise(vap_n = n(), .groups = "drop")
episode_summary_vap <- vap_summary %>%
  full_join(summarise_pathogen(vap_episode_icu), by=c("siteid","country_income")) %>%
  mutate(across(all_of(pathogen_list), ~replace_na(., 0L)))

bsi_summary_all <- bsi_episode_all %>%
  group_by(siteid, country_income) %>%
  summarise(bsi_n = n(), .groups = "drop")
episode_summary_bsi_all <- bsi_summary_all %>%
  full_join(summarise_pathogen(bsi_episode_all), by=c("siteid","country_income")) %>%
  mutate(across(all_of(pathogen_list), ~replace_na(., 0L)))

bsi_summary_icu <- bsi_episode_icu %>%
  group_by(siteid, country_income) %>%
  summarise(bsi_n = n(), .groups = "drop")
episode_summary_bsi_icu <- bsi_summary_icu %>%
  full_join(summarise_pathogen(bsi_episode_icu), by=c("siteid","country_income")) %>%
  mutate(across(all_of(pathogen_list), ~replace_na(., 0L)))

# Denominator builder
prepare_bed_period <- function(df, study_period_updated,
                               icu_multiplier = 5,
                               imcu_multiplier = 2,
                               cap_other_beds = 100,
                               cap_other = TRUE,
                               icu_imcu_only = FALSE) {
  df_bed <- df %>%
    select(Site, Ward, Beds, ward_cal) %>%
    distinct() %>%
    filter(!if_all(everything(), ~ is.na(.) | . == "")) %>%
    mutate(Beds = as.numeric(str_replace_all(Beds, "[^0-9.]", "")))
  
  if (icu_imcu_only) {
    df_bed <- df_bed %>% filter(ward_cal %in% c("ICU","IMCU"))
  } else {
    df_bed <- df_bed %>%
      mutate(Beds = case_when(
        ward_cal == "ICU"  ~ Beds * icu_multiplier,
        ward_cal == "IMCU" ~ Beds * imcu_multiplier,
        TRUE ~ ifelse(cap_other & ward_cal == "Other" & Beds > cap_other_beds, cap_other_beds, Beds)
      ))
  }
  
  # fill missing beds by site mean
  df_bed <- df_bed %>%
    group_by(Site) %>%
    mutate(Beds = ifelse(is.na(Beds), round(mean(Beds, na.rm = TRUE)), Beds)) %>%
    ungroup()
  
  period_all <- study_period_updated %>%
    transmute(Site, period = period_valid_final) %>%
    left_join(df_bed, by = "Site") %>%
    mutate(period_bed = ifelse(is.na(period) | is.na(Beds), 0, period * Beds)) %>%
    group_by(Site) %>%
    summarise(total_period_bed = sum(period_bed, na.rm = TRUE), .groups = "drop") %>%
    filter(total_period_bed > 0)
  
  list(period_all = period_all)
}

# VAP: ICU+IMCU original beds
res_vap     <- prepare_bed_period(df, study_period_updated, icu_imcu_only = TRUE)
# BSI-ALL: ICU×5, IMCU×2, Other (capped)
res_bsi_all <- prepare_bed_period(df, study_period_updated, icu_imcu_only = FALSE)
# BSI-ICU: ICU+IMCU original beds
res_bsi_icu <- prepare_bed_period(df, study_period_updated, icu_imcu_only = TRUE)


# Merge denominators + events
normalize_meta_df <- function(d) {
  d %>%
    left_join(site_map, by="siteid") %>%
    mutate(
      site_label = factor(site_label, levels = site_map$site_label),
      country_income = factor(country_income, levels = c("High income","Upper middle income","Lower middle income"))
    ) %>%
    arrange(siteid)
}

# Totals
df_vap      <- res_vap$period_all     %>% rename(siteid = Site) %>%
  left_join(episode_summary_vap, by = "siteid")     %>% normalize_meta_df()

df_bsi_all  <- res_bsi_all$period_all %>% rename(siteid = Site) %>%
  left_join(episode_summary_bsi_all, by = "siteid") %>% normalize_meta_df()

df_bsi_icu  <- res_bsi_icu$period_all %>% rename(siteid = Site) %>%
  left_join(episode_summary_bsi_icu, by = "siteid") %>% normalize_meta_df()

# Pathogens
df_vap_patho <- df_vap
df_bsi_patho_all  <- df_bsi_all
df_bsi_patho_icu  <- df_bsi_icu

# ----------------------------------------------------------------------
# Forest plot 
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_open_pdf <- function(file, width = 11.5, height = 14.5, family = "Times") {
  ok <- FALSE
  if (requireNamespace("Cairo", quietly = TRUE)) {
    try({ Cairo::CairoPDF(file, width = width, height = height, family = family); ok <- TRUE }, silent = TRUE)
  }
  if (!ok) grDevices::pdf(file, width = width, height = height, family = family)
}

fit_group_with_weights <- function(dat, event_var, time_var, scale = 10000) {
  if (nrow(dat) == 0) return(NULL)
  m <- meta::metarate(
    event = dat[[event_var]],
    time  = dat[[time_var]],
    studlab = dat$site_label,
    sm = "IRLN",
    method.tau = "REML",
    random = TRUE,
    method.random.ci = "HK",
    overall = TRUE
  )
  grp <- tibble(
    est = exp(m$TE.random) * scale,
    lcl = exp(m$lower.random) * scale,
    ucl = exp(m$upper.random) * scale,
    I2 = m$I2, tau2 = m$tau2, Q = m$Q, pQ = m$pval.Q
  )
  w_raw <- m$w.random
  s <- sum(w_raw, na.rm = TRUE)
  w_pct <- if (is.finite(s) && s > 0) 100 * w_raw / s else rep(NA_real_, length(w_raw))
  w <- tibble(studlab = m$studlab, weight = as.numeric(w_pct))
  list(grp = grp, weights = w)
}

build_forest_rows_rich <- function(data, event_var,
                                   time_var = "total_period_bed",
                                   lab_var = "site_label",
                                   group_var  = "country_income",
                                   scale = 10000,
                                   order_desc = TRUE) {
  df <- data %>%
    mutate(.event = .data[[event_var]],
           .time  = .data[[time_var]],
           .lab   = .data[[lab_var]],
           .grp   = .data[[group_var]]) %>%
    select(.grp, .lab, .event, .time, everything())
  
  grp_levels <- if (is.factor(df$.grp)) levels(df$.grp) else sort(unique(df$.grp))
  rows_list <- list()
  
  extra_gap_for <- function(lbl) {
    if (lbl %in% c("Upper middle income", "Lower middle income")) return(0.8)
    if (lbl %in% c("Overall")) return(1.0)
    0
  }
  
  first_group <- TRUE 
  
  for (g in grp_levels) {
    dfg <- df %>% filter(.grp == g)
    if (nrow(dfg) == 0) next
    
    fg <- fit_group_with_weights(dfg, event_var, time_var, scale)
    grp_stat <- fg$grp %>% mutate(.grp = g)
    
    l <- ifelse(dfg$.event == 0, 0, stats::qchisq(0.025, 2*dfg$.event) / (2*dfg$.time))
    u <- stats::qchisq(0.975, 2*(dfg$.event + 1)) / (2*dfg$.time)
    r <- dfg$.event / dfg$.time
    rc <- tibble(rate = r*scale, lcl = l*scale, ucl = u*scale)
    
    dfg <- dfg %>%
      left_join(fg$weights, by = c(".lab" = "studlab")) %>%
      bind_cols(rc) %>%
      arrange(if (order_desc) desc(rate) else rate)
    
    
    gap <- extra_gap_for(as.character(g))
    if (!first_group && gap > 0) {
      rows_list[[length(rows_list)+1]] <- tibble(
        type = "spacer", label = "", .grp = g,
        rate = NA_real_, lcl = NA_real_, ucl = NA_real_,
        events = NA_integer_, time = NA_real_, weight = NA_real_,
        y_offset = gap
      )
    }
    
    rows_list[[length(rows_list)+1]] <- tibble(
      type = "header", label = as.character(g), .grp = g,
      rate = NA_real_, lcl = NA_real_, ucl = NA_real_,
      events = NA_integer_, time = NA_real_, weight = NA_real_,
      y_offset = 1
    )
    
    rows_list[[length(rows_list)+1]] <- dfg %>%
      transmute(
        type = "study", label = .lab, .grp = .grp,
        rate, lcl, ucl, events = .event, time = .time, weight = weight,
        y_offset = 1
      )
    
    rows_list[[length(rows_list)+1]] <- tibble(
      type = "pooled_group",
      label = sprintf(
        "\nRandom-effects (REML)\n(I\u00B2=%.1f%%, \u03C4\u00B2=%.3f, Q=%.2f, p=%s)",
        grp_stat$I2, grp_stat$tau2, grp_stat$Q,
        ifelse(is.na(grp_stat$pQ), "NA", formatC(grp_stat$pQ, format="f", digits=4))
      ),
      .grp = g,
      rate = grp_stat$est, lcl = grp_stat$lcl, ucl = grp_stat$ucl,
      events = NA_integer_, time = NA_real_, weight = NA_real_,
      y_offset = 1
    )
    
    first_group <- FALSE  
  }
  
  
  sg_all <- meta::metarate(
    event = df$.event, time = df$.time, studlab = df$.lab,
    sm = "IRLN", method.tau = "REML", random = TRUE,
    method.random.ci = "HK", overall = TRUE
  )
  overall_est <- tibble(
    est = exp(sg_all$TE.random) * scale,
    lcl = exp(sg_all$lower.random) * scale,
    ucl = exp(sg_all$upper.random) * scale,
    I2 = sg_all$I2, tau2 = sg_all$tau2, Q = sg_all$Q, pQ = sg_all$pval.Q
  )
  
  gap_overall <- extra_gap_for("Overall")
  if (gap_overall > 0) {
    rows_list[[length(rows_list)+1]] <- tibble(
      type = "spacer", label = "", .grp = "Overall",
      rate = NA_real_, lcl = NA_real_, ucl = NA_real_,
      events = NA_integer_, time = NA_real_, weight = NA_real_,
      y_offset = gap_overall
    )
  }
  
  rows_list[[length(rows_list)+1]] <- tibble(
    type = "header", label = "Overall", .grp = "Overall",
    rate = NA_real_, lcl = NA_real_, ucl = NA_real_,
    events = NA_integer_, time = NA_real_, weight = NA_real_,
    y_offset = 1
  )
  
  rows_list[[length(rows_list)+1]] <- tibble(
    type = "pooled_overall",
    label = sprintf(
      "Random-effects (REML)\n(I\u00B2=%.1f%%, \u03C4\u00B2=%.3f, Q=%.2f, p=%g)",
      overall_est$I2, overall_est$tau2, overall_est$Q, overall_est$pQ
    ),
    .grp = "Overall",
    rate = overall_est$est, lcl = overall_est$lcl, ucl = overall_est$ucl,
    events = NA_integer_, time = NA_real_, weight = NA_real_,
    y_offset = 1
  )
  
  rows <- bind_rows(rows_list)
  rows$y <- rev(cumsum(rows$y_offset))
  list(rows = rows)
}


plot_forest_ir_rich <- function(layout,
                                xlab = "Incidence (per 10000 patient-days)",
                                time_label = "Time (bed-days)",
                                xlim_upper = 50,
                                diamond_h = 0.35,
                                font_study = "Times",
                                font_num   = "Times") {
  rows <- layout$rows
  
  fmt_num <- function(x, d=2) ifelse(is.na(x), "", formatC(x, digits=d, format="f"))
  right_col <- dplyr::case_when(
    rows$type %in% c("study","pooled_group","pooled_overall") ~
      paste0(fmt_num(rows$rate)," [",fmt_num(rows$lcl),", ",fmt_num(rows$ucl),"]"),
    TRUE ~ ""
  )
  left_events <- ifelse(is.na(rows$events), "", as.character(rows$events))
  left_time   <- ifelse(is.na(rows$time),   "", format(round(rows$time), big.mark=""))
  
  diamonds <- rows %>%
    dplyr::filter(type %in% c("pooled_group","pooled_overall")) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::rowwise() %>%
    do({
      r <- .
      tibble::tibble(
        id = r$id,
        x = c(r$lcl, r$rate, r$ucl, r$rate),
        y = c(r$y,   r$y + diamond_h, r$y, r$y - diamond_h),
        label = r$label
      )
    }) %>% dplyr::bind_rows()
  
  xmin <- -0.80 * xlim_upper
  xmax <-  1.10 * xlim_upper
  x_left_study  <- -0.72 * xlim_upper
  x_left_events <- -0.30 * xlim_upper
  x_left_time   <- -0.06 * xlim_upper
  x_right_rate  <-  1.06 * xlim_upper
  
  y_max  <- max(rows$y)
  y_head <- y_max + 1.3
  
  axis_breaks <- pretty(c(0, xlim_upper))
  axis_breaks <- axis_breaks[axis_breaks >= 0]
  
  ggplot() +
    geom_errorbarh(data = dplyr::filter(rows, type == "study"),
                   aes(y = y, xmin = pmax(lcl, 0), xmax = pmax(ucl, 0)),
                   height = 0.25, linewidth = 0.3) +
    geom_point(
      data = dplyr::filter(rows, type == "study"),
      aes(x = pmax(rate, 0), y = y),
      shape = 16, size = 0.3
    ) +
    geom_polygon(data = diamonds %>% dplyr::mutate(x = pmax(x, 0)),
                 aes(x = x, y = y, group = id), alpha = 0.5) +
    geom_text(data = rows,
              aes(x = x_left_study, y = y, label = label),
              hjust = 0, size = 3.4, family = font_study,
              lineheight = 1.08, na.rm = TRUE, fontface = "bold") +
    geom_text(
      data = rows, aes(x = x_left_events, y = y), label = left_events,
      hjust = 1, size = 3.4, family = font_num, na.rm = TRUE, fontface = "bold"
    ) +
    geom_text(
      data = rows, aes(x = x_left_time, y = y), label = left_time,
      hjust = 1, size = 3.4, family = font_num, na.rm = TRUE, fontface = "bold"
    ) +
    geom_text(
      data = rows, aes(x = x_right_rate, y = y), label = right_col,
      hjust = 0, size = 3.4, family = font_num, na.rm = TRUE, fontface = "bold"
    ) +
    annotate("text", x = x_left_study,  y = y_head, label = "Study",
             hjust = 0, size = 3.4, family = font_study, fontface = "bold") +
    annotate("text", x = x_left_events, y = y_head, label = "Events",
             hjust = 1, size = 3.4, family = font_study, fontface = "bold") +
    annotate("text", x = x_left_time,   y = y_head, label = time_label,
             hjust = 1, size = 3.4, family = font_study, fontface = "bold") +
    annotate("text", x = x_right_rate,  y = y_head, label = "Incidence [95% CI]",
             hjust = 0, size = 3.4, family = font_study, fontface = "bold") +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0, y_head + 1), expand = c(0, 0)) +
    scale_x_continuous(
      xlab, limits = c(xmin, xmax),
      breaks = axis_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    geom_segment(aes(x = -0.05, xend = xlim_upper+0.05, y = 0, yend = 0),
                 inherit.aes = FALSE, linewidth = 0.3) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length.x = grid::unit(3, "pt"),
      axis.title.x = element_text(colour = "black", size = 9,
                                  family = "Times", face = "bold",
                                  margin = margin(t = 5), hjust = 0.7),
      axis.text.x = element_text(colour = "black", size = 9,
                                 family = "Times", face = "bold"),
      plot.margin = margin(10, 180, 5, 20)
    )
}

# Drop entire income groups if pooled upper CI > cutoff
filter_groups_by_upper_ci <- function(data,
                                      event_var,
                                      time_var   = "total_period_bed",
                                      group_var  = "country_income",
                                      scale      = 10000,
                                      upper_cut  = 200) {
  d <- data %>% dplyr::filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0)
  if (nrow(d) == 0) return(d[0, ])
  
  gm <- d %>%
    dplyr::filter(!is.na(.data[[group_var]]), .data[[group_var]] != "") %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::group_map(~{
      xx <- .x
      u  <- Inf
      if (nrow(xx) > 0) {
        m <- tryCatch(
          meta::metarate(
            event   = xx[[event_var]],
            time    = xx[[time_var]],
            studlab = xx$site_label,
            sm = "IRLN",
            method.tau = "REML",
            random = TRUE,
            method.random.ci = "HK",
            overall = TRUE
          ),
          error = function(e) NULL
        )
        if (!is.null(m)) u <- exp(m$upper.random) * scale
      }
      tibble::tibble(grp = .y[[group_var]][1], ucl = u)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::ungroup()
  
  keep <- gm %>% dplyr::filter(ucl <= upper_cut) %>% dplyr::pull(grp)
  dropped <- setdiff(unique(d[[group_var]]), keep)
  if (length(dropped) > 0) {
    message("Dropped income groups (pooled upper CI > ", upper_cut, "): ",
            paste(dropped, collapse = ", "))
  }
  d %>% dplyr::filter(.data[[group_var]] %in% keep)
}


run_meta_analysis <- function(data,
                              event_var,
                              time_var = "total_period_bed",
                              scale = 10000,
                              xlab = "Incidence (per 10000 patient-days)",
                              time_label = "Time (bed-days)",
                              xlim_upper = 50,
                              pdf_file,
                              width = 11.5,
                              height = 14.5,
                              font_family = "Times",
                              upper_ci_cutoff = 1000000) {
  d <- data %>% dplyr::filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0)
  if (nrow(d) == 0) {
    message("No data to analyze for: ", event_var, " (", pdf_file, ")")
    return(invisible(NULL))
  }
  
  # Drop income levels whose pooled upper CI exceeds cutoff
  d <- filter_groups_by_upper_ci(
    data      = d,
    event_var = event_var,
    time_var  = time_var,
    group_var = "country_income",
    scale     = scale,
    upper_cut = upper_ci_cutoff
  )
  if (nrow(d) == 0) {
    message("All groups removed by pooled upper CI filter (> ", upper_ci_cutoff, "): ", pdf_file)
    return(invisible(NULL))
  }
  
  layout <- build_forest_rows_rich(
    data = d,
    event_var= event_var,
    time_var = time_var,
    lab_var  = "site_label",
    group_var= "country_income",
    scale = scale
  )
  safe_open_pdf(pdf_file, width = width, height = height, family = font_family)
  p <- plot_forest_ir_rich(layout,
                           xlab = xlab,
                           time_label = time_label,
                           xlim_upper = xlim_upper,
                           font_study = font_family,
                           font_num = font_family)
  grid::grid.newpage(); grid::grid.draw(ggplotGrob(p))
  grDevices::dev.off()
}


# BSI ALL wards
run_meta_analysis(
  data        = df_bsi_all %>% filter(bsi_n > 0),
  event_var   = "bsi_n",
  xlab        = "BSI incidence (per 10000 patient-days)",
  xlim_upper  = 50,
  pdf_file    = "output/figure/incidence_meta_BSI_ALL.pdf",
  upper_ci_cutoff = 1000000
)

# BSI ICU-only
run_meta_analysis(
  data        = df_bsi_icu %>% filter(bsi_n > 0),
  event_var   = "bsi_n",
  xlab        = "BSI incidence (per 10000 ICU-days)",
  xlim_upper  = 600,
  pdf_file    = "output/figure/incidence_meta_BSI_ICU.pdf",
  upper_ci_cutoff = 1000000
)

# VAP
run_meta_analysis(
  data        = df_vap %>% filter(vap_n > 0),
  event_var   = "vap_n",
  xlab        = "VAP incidence (per 10000 ICU-days)",
  xlim_upper  = 600,
  pdf_file    = "output/figure/incidence_meta_VAP.pdf",
  upper_ci_cutoff = 1000000
)

# VAP pathogens (ICU-only denominator)
pathogen_xlim_vap <- c(CRA=300, `3GCRE`=150, CRE=200, CRP=70, VRE=20, MRSA=120)
pathogen_height_vap <- c(CRA=13, `3GCRE`=13, CRE=13, CRP=11, VRE=3, MRSA=8) 

for (var in names(pathogen_xlim_vap)) {
  df_tmp <- df_vap_patho %>% filter(.data[[var]] > 0)
  run_meta_analysis(
    data        = df_tmp,
    event_var   = var,
    xlab        = "VAP incidence (per 10000 ICU-days)",
    xlim_upper  = pathogen_xlim_vap[[var]],
    pdf_file    = paste0("output/figure/incidence_meta_VAP_", var, ".pdf"),
    upper_ci_cutoff = 1000000,
    height = pathogen_height_vap[[var]]
  )
}

# BSI pathogens - ALL wards denominator
pathogen_xlim_bsi <- c(CRA=20, `3GCRE`=15, CRE=25, CRP=3, VRE=5, MRSA=6)
pathogen_height_bsi <- c(CRA=14.5, `3GCRE`=14.5, CRE=13, CRP=11, VRE=8, MRSA=12) 

for (var in names(pathogen_xlim_bsi)) {
  df_tmp <- df_bsi_patho_all %>% filter(.data[[var]] > 0)
  run_meta_analysis(
    data        = df_tmp,
    event_var   = var,
    xlab        = "BSI incidence (per 10000 patient-days)",
    xlim_upper  = pathogen_xlim_bsi[[var]],
    pdf_file    = paste0("output/figure/incidence_meta_BSI_ALL_", var, ".pdf"),
    upper_ci_cutoff = 1000000,
    height = pathogen_height_bsi[[var]]
  )
}

# BSI pathogens - ICU-only denominator
pathogen_xlim_bsi_ICU <- c(CRA=150, `3GCRE`=250, CRE=250, CRP=30, VRE=10, MRSA=50)
pathogen_height_bsi_ICU <- c(CRA=13, `3GCRE`=13, CRE=13, CRP=11, VRE=7, MRSA=6.5) 

for (var in names(pathogen_xlim_bsi_ICU)) {
  df_tmp <- df_bsi_patho_icu %>% filter(.data[[var]] > 0)
  run_meta_analysis(
    data        = df_tmp,
    event_var   = var,
    xlab        = "BSI incidence (per 10000 ICU-days)",
    xlim_upper  = pathogen_xlim_bsi_ICU[[var]],
    pdf_file    = paste0("output/figure/incidence_meta_BSI_ICU_", var, ".pdf"),
    upper_ci_cutoff = 1000000,
    height = pathogen_height_bsi_ICU[[var]]
  )
}


# ---------------------------------------------------
# Save outcomes
.meta_one <- function(d, event_var, time_var = "total_period_bed", scale = 10000) {
  d <- d %>% filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0)
  if (nrow(d) == 0) return(tibble(
    k = 0, events = NA_integer_, time = NA_real_,
    est = NA_real_, lcl = NA_real_, ucl = NA_real_,
    I2 = NA_real_, tau2 = NA_real_, Q = NA_real_, pQ = NA_real_
  ))
  m <- meta::metarate(
    event = d[[event_var]],
    time = d[[time_var]],
    studlab = d$site_label,
    sm = "IRLN",
    method.tau = "REML",
    random = TRUE,
    method.random.ci = "HK",
    overall = TRUE
  )
  tibble(
    k = length(m$studlab),
    events = sum(d[[event_var]], na.rm = TRUE),
    time = sum(d[[time_var]], na.rm = TRUE),
    est = exp(m$TE.random) * scale,
    lcl = exp(m$lower.random) * scale,
    ucl = exp(m$upper.random) * scale,
    I2 = m$I2,
    tau2 = m$tau2,
    Q = m$Q,
    pQ = m$pval.Q
  )
}

compute_meta_tables <- function(data, event_var,
                                group_var = "country_income",
                                time_var = "total_period_bed",
                                scale = 10000,
                                upper_ci_cutoff = 1000000) {
  # Apply the same group-removal rule to the table inputs
  d0 <- data %>%
    dplyr::filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0) %>%
    filter_groups_by_upper_ci(
      event_var = event_var,
      time_var  = time_var,
      group_var = group_var,
      scale     = scale,
      upper_cut = upper_ci_cutoff
    ) %>%
    dplyr::filter(!is.na(.data[[group_var]]) & .data[[group_var]] != "") %>%
    dplyr::mutate(.grp = .data[[group_var]])
  
  if (nrow(d0) == 0) {
    return(list(
      table = tibble(Group=character(), k=integer(), events=integer(), time=double(),
                     est=double(), lcl=double(), ucl=double(), rate_ci=character(),
                     I2=double(), tau2=double(), Q=double(), pQ=double()),
      Qb = NA_real_, pQb = NA_real_
    ))
  }
  
  per_grp <- d0 %>%
    dplyr::group_by(.grp) %>%
    dplyr::group_modify(~ .meta_one(.x, event_var, time_var, scale)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(Group = .grp) %>%
    dplyr::mutate(Group = as.character(Group))
  
  overall_row <- .meta_one(d0, event_var, time_var, scale) %>%
    dplyr::mutate(Group = "Overall") %>%
    dplyr::select(Group, dplyr::everything())
  
  table_out <- dplyr::bind_rows(per_grp, overall_row) %>%
    dplyr::mutate(rate_ci = dplyr::case_when(
      is.na(est) ~ NA_character_,
      TRUE ~ sprintf("%.2f [%.2f, %.2f]", est, lcl, ucl)
    )) %>%
    dplyr::select(Group, k, events, time, est, lcl, ucl, rate_ci, I2, tau2, Q, pQ)
  
  m_sub <- meta::metarate(
    event   = d0[[event_var]],
    time    = d0[[time_var]],
    studlab = d0$site_label,
    subgroup = d0[[group_var]],
    sm = "IRLN",
    method.tau = "REML",
    random = TRUE,
    method.random.ci = "HK",
    overall = TRUE
  )
  Qb <- suppressWarnings(tryCatch(m_sub$Q.b.random, error = function(e) NA_real_))
  pQb <- suppressWarnings(tryCatch(m_sub$pval.Q.b,   error = function(e) NA_real_))
  
  list(table = table_out, Qb = Qb, pQb = pQb)
}

# Tables: totals
res_bsi_all_tab <- compute_meta_tables(df_bsi_all %>% filter(bsi_n > 0), "bsi_n")
res_bsi_icu_tab <- compute_meta_tables(df_bsi_icu %>% filter(bsi_n > 0), "bsi_n")
res_vap_tab <- compute_meta_tables(df_vap %>% filter(vap_n > 0), "vap_n")

# Tables: pathogens (VAP / BSI-ALL / BSI-ICU)
res_vap_patho_tabs <- purrr::map(pathogen_list,
                                 ~ compute_meta_tables(df_vap_patho %>% filter(.data[[.x]] > 0), .x))
names(res_vap_patho_tabs) <- pathogen_list

res_bsi_all_patho_tabs <- purrr::map(pathogen_list,
                                     ~ compute_meta_tables(df_bsi_patho_all %>% filter(.data[[.x]] > 0), .x))
names(res_bsi_all_patho_tabs) <- pathogen_list

res_bsi_icu_patho_tabs <- purrr::map(pathogen_list,
                                     ~ compute_meta_tables(df_bsi_patho_icu %>% filter(.data[[.x]] > 0), .x))
names(res_bsi_icu_patho_tabs) <- pathogen_list

# Write Excel
wb <- openxlsx::createWorkbook()

.add_sheet <- function(wb, sheet_name, tbl) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = tbl)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(tbl), widths = "auto")
}

# totals
.add_sheet(wb, "BSI_ALL", res_bsi_all_tab$table)
.add_sheet(wb, "BSI_ICU", res_bsi_icu_tab$table)
.add_sheet(wb, "VAP", res_vap_tab$table)

# pathogens — VAP
for (nm in names(res_vap_patho_tabs)) .add_sheet(wb, paste0("VAP_", nm), res_vap_patho_tabs[[nm]]$table)
# pathogens — BSI ALL
for (nm in names(res_bsi_all_patho_tabs)) .add_sheet(wb, paste0("BSI_ALL_", nm), res_bsi_all_patho_tabs[[nm]]$table)
# pathogens — BSI ICU
for (nm in names(res_bsi_icu_patho_tabs)) .add_sheet(wb, paste0("BSI_ICU_", nm), res_bsi_icu_patho_tabs[[nm]]$table)

subgroup_tests <- dplyr::bind_rows(
  tibble::tibble(Analysis = "BSI_ALL", Qb = res_bsi_all_tab$Qb, pQb = res_bsi_all_tab$pQb),
  tibble::tibble(Analysis = "BSI_ICU", Qb = res_bsi_icu_tab$Qb, pQb = res_bsi_icu_tab$pQb),
  tibble::tibble(Analysis = "VAP", Qb = res_vap_tab$Qb, pQb = res_vap_tab$pQb),
  purrr::imap_dfr(res_vap_patho_tabs, ~ tibble::tibble(Analysis = paste0("VAP_", .y), Qb = .x$Qb, pQb = .x$pQb)),
  purrr::imap_dfr(res_bsi_all_patho_tabs, ~ tibble::tibble(Analysis = paste0("BSI_ALL_", .y), Qb = .x$Qb, pQb = .x$pQb)),
  purrr::imap_dfr(res_bsi_icu_patho_tabs, ~ tibble::tibble(Analysis = paste0("BSI_ICU_", .y), Qb = .x$Qb, pQb = .x$pQb))
)
.add_sheet(wb, "Subgroup_tests", subgroup_tests)

xlsx_file <- "output/table/incidence_meta_summaries.xlsx"
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)


