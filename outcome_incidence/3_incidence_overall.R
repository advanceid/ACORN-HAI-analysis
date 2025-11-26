# Clear workspace
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    readxl, dplyr, tidyr, ggplot2, forcats,
    Cairo, ggh4x, ggpp, purrr, scales, stringr
  )
})

# Inputs / outputs
xlsx_in  <- "output/incidence_meta_summaries.xlsx"
pdf_out  <- "pdf/summary/incidence_summary.pdf"

# Figure setup
base_family <- "Times"
fig_width   <- 7
fig_height  <- 4

# Only three pathogens
categories <- c("CRA", "3GCRE", "CRE")

analysis_map <- tibble::tribble(
  ~sheet_key, ~analysis_label,
  "VAP",      "VAP\n(per 10000 ICU-days)",
  "BSI_ICU",  "Hospital-acquired BSI\n(per 10000 ICU-days)",
  "BSI_ALL",  "Hospital-acquired BSI\n(per 10000 patient-days)"
)

pathogen_sheet_map <- tibble::tribble(
  ~sheet_key, ~sheet_name,
  "VAP",      "VAP_pathogens_all",
  "BSI_ICU",  "BSI_ICU_pathogens_all",
  "BSI_ALL",  "BSI_ALL_pathogens_all"
)

# Group order for y-axis
group_levels <- c("High income", "Upper middle income", "Lower middle income", "Overall")

# Helpers
.read_sheet <- function(sheet_name) {
  dat <- suppressWarnings(readxl::read_excel(xlsx_in, sheet = sheet_name))
  names(dat) <- tolower(names(dat))
  dat
}

.sheet_exists <- function(sheet_name) {
  sheet_name %in% readxl::excel_sheets(xlsx_in)
}

# Read one block (one analysis × one pathogen)
.fetch_block <- function(sheet_key, analysis_label, category) {
  preferred <- pathogen_sheet_map %>%
    dplyr::filter(sheet_key == !!sheet_key) %>%
    dplyr::pull(sheet_name)
  stopifnot(length(preferred) == 1)
  
  if (.sheet_exists(preferred)) {
    dat <- .read_sheet(preferred)
    req_cols <- c("pathogen", "group", "est", "lcl", "ucl")
    if (!all(req_cols %in% names(dat))) {
      stop("Consolidated sheet '", preferred, "' is missing required columns: ",
           paste(setdiff(req_cols, names(dat)), collapse = ", "))
    }
    out <- dat %>%
      dplyr::filter(.data[["pathogen"]] == category) %>%
      dplyr::mutate(
        analysis = analysis_label,
        category = category,
        group    = factor(as.character(group), levels = group_levels)
      ) %>%
      dplyr::select(analysis, category, group, est, lcl, ucl)
    return(out)
  } else {
    # Fallback to sheet per pathogen, e.g. "VAP_CRA"
    fallback_sheet <- paste0(sheet_key, "_", category)
    if (!.sheet_exists(fallback_sheet)) {
      stop("Neither consolidated sheet '", preferred,
           "' nor fallback sheet '", fallback_sheet, "' exists in: ", xlsx_in)
    }
    dat <- .read_sheet(fallback_sheet)
    req_cols <- c("group", "est", "lcl", "ucl")
    if (!all(req_cols %in% names(dat))) {
      stop("Fallback sheet '", fallback_sheet, "' is missing required columns: ",
           paste(setdiff(req_cols, names(dat)), collapse = ", "))
    }
    out <- dat %>%
      dplyr::mutate(
        analysis = analysis_label,
        category = category,
        group    = factor(as.character(group), levels = group_levels)
      ) %>%
      dplyr::select(analysis, category, group, est, lcl, ucl)
    return(out)
  }
}

# Plot data
df_plot <- purrr::map_dfr(
  seq_len(nrow(analysis_map)),
  function(i) {
    ak <- analysis_map$sheet_key[i]
    al <- analysis_map$analysis_label[i]
    purrr::map_dfr(categories, ~ .fetch_block(ak, al, .x))
  }
)

# Ensure numeric
df_plot <- df_plot %>%
  mutate(across(c(est, lcl, ucl), as.numeric))

# Facet order
df_plot <- df_plot %>%
  mutate(
    category = factor(category, levels = categories),
    analysis = factor(analysis, levels = analysis_map$analysis_label)
  )

# y order (top to bottom)
lvl_y <- rev(c("High income","Upper middle income","Lower middle income","Overall"))

# Per-panel x upper limits
df_all <- df_plot %>%
  mutate(
    group   = factor(as.character(group), levels = lvl_y),
    y_num   = as.numeric(group),
    est_p   = est,
    lcl_p   = lcl,
    ucl_p   = ucl,
    x_cap   = dplyr::case_when(
      analysis == analysis_map$analysis_label[3] ~ 4, 
      TRUE ~ 20
    )
  )

df_income <- df_all %>% filter(group != "Overall")
df_over   <- df_all %>% filter(group == "Overall")

# CI segment endpoints + truncation flag
df_income <- df_income %>%
  mutate(
    need_break = ucl_p > x_cap,
    xmin_plot  = pmax(lcl_p, 0),
    xmax_plot  = pmin(ucl_p, x_cap),  
    x_point    = pmin(est_p, x_cap) 
  )

df_over <- df_over %>%
  mutate(
    need_break = ucl_p > x_cap,
    xmin_plot  = pmax(lcl_p, 0),
    xmax_plot  = pmin(ucl_p, x_cap),
    x_point    = pmin(est_p, x_cap)
  )

cap_half <- 0.18

# Dummy data to enforce panel-specific x limits
df_limits <- tidyr::expand_grid(
  analysis = factor(analysis_map$analysis_label, levels = analysis_map$analysis_label),
  category = factor(categories, levels = categories),
  y_num    = 1,
  which    = c("min", "max")
) %>%
  mutate(
    x = dplyr::case_when(
      which == "min" ~ 0,
      which == "max" & analysis == analysis_map$analysis_label[3] ~ 4,
      which == "max" ~ 20
    )
  )

# Plot
p <- ggplot() +
# -------------------------
# Income groups (non-Overall)
# -------------------------
geom_segment(
  data = df_income,
  aes(x = xmin_plot, xend = xmax_plot, y = y_num, yend = y_num),
  linewidth = 0.25, colour = "black"
) +
  # Left caps
  geom_segment(
    data = df_income,
    aes(x = xmin_plot, xend = xmin_plot,
        y = y_num - cap_half, yend = y_num + cap_half),
    linewidth = 0.25, colour = "black"
  ) +
  # Right caps for non-truncated CIs
  geom_segment(
    data = df_income %>% filter(!need_break),
    aes(x = xmax_plot, xend = xmax_plot,
        y = y_num - cap_half, yend = y_num + cap_half),
    linewidth = 0.25, colour = "black"
  ) +
  geom_text(
    data = df_income %>% 
      filter(need_break, analysis != analysis_map$analysis_label[3]),
    aes(x = xmax_plot - 0.15, y = y_num),
    label  = "//",
    family = base_family,
    fontface = "bold",
    size   = 2.6,
    hjust  = 0,
    vjust  = 0.4,
    colour = "black"
  ) +
  geom_point(
    data = df_income %>%
      filter(need_break, analysis == analysis_map$analysis_label[3]),
    aes(x = xmax_plot+0.05, y = y_num),
    shape = 2,    
    size  = 1.5,
    stroke = 0.5,
    colour = "black"
  ) +
  geom_point(
    data = df_income,
    aes(x = x_point, y = y_num),
    shape = 19, size = 0.8, colour = "black"
  ) +
  
# -------------------------
# Overall group
# -------------------------
geom_segment(
  data = df_over,
  aes(x = xmin_plot, xend = xmax_plot, y = y_num, yend = y_num),
  linewidth = 0.25, colour = "black"
) +
  # Left caps
  geom_segment(
    data = df_over,
    aes(x = xmin_plot, xend = xmin_plot,
        y = y_num - cap_half, yend = y_num + cap_half),
    linewidth = 0.25, colour = "black"
  ) +
  # Right caps for non-truncated CIs
  geom_segment(
    data = df_over %>% filter(!need_break),
    aes(x = xmax_plot, xend = xmax_plot,
        y = y_num - cap_half, yend = y_num + cap_half),
    linewidth = 0.25, colour = "black"
  ) +
  geom_text(
    data = df_over %>% 
      filter(need_break, analysis != analysis_map$analysis_label[3]),
    aes(x = xmax_plot - 0.15, y = y_num),
    label  = "//",
    family = base_family,
    fontface = "bold",
    size   = 2.6,
    hjust  = 0,
    vjust  = 0.4,
    colour = "black"
  ) +
  geom_point(
    data = df_over %>% 
      filter(need_break, analysis == analysis_map$analysis_label[3]),
    aes(x = xmax_plot+0.05, y = y_num),
    shape = 2,
    size  = 1.5,
    stroke = 0.5,
    colour = "black"
  ) +
  geom_point(
    data = df_over,
    aes(x = x_point, y = y_num),
    shape = 18, size = 2, colour = "black"
  ) +
  
  # Dummy points to enforce x-range per facet
  geom_blank(
    data = df_limits,
    aes(x = x, y = y_num)
  ) +
  # y-axis: numeric -> income labels
  scale_y_continuous(
    name   = NULL,
    breaks = seq_along(lvl_y),
    labels = lvl_y,
    limits = c(0.5, length(lvl_y) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_minimal(base_family = base_family, base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y  = element_text(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey80", color = NA),
    legend.position = "none",
    plot.margin = margin(6,8,6,8),
    panel.spacing.x = unit(1, "lines")
  ) +
  facet_grid2(vars(category), vars(analysis), scales = "free_x", independent = "x") 

p

# Save PDF
Cairo::CairoPDF(
  file   = pdf_out,
  width  = fig_width,
  height = fig_height,
  family = base_family,
  onefile = TRUE
)
print(p)
dev.off()
