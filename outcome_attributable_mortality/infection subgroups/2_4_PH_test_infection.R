# Clear workspace
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, survival, gt, gtsummary, purrr, extrafont)
})

# Set working directory and load fonts
setwd("./")
loadfonts()

# Load data
df_all <- readRDS("data/att_first28_infection.RData")

# Define pathogen labels & output paths
bug_labels <- c(
  "Acinetobacter spp.",
  "Enterobacterales (3GCRE)",
  "Enterobacterales (CRE)",
  "Pseudomonas spp.",
  "Enterococcus spp.",
  "Staphylococcus aureus",
  "MDR GNB"
)
resist_vars  <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth", "mdr_gnb")
output_paths <- paste0(
  "output/table/cox_zph_att28_infection_",
  c("01_cra","02_3gcre","03_cre","04_crp","05_vre","06_mrsa", "07_mdr"),
  ".html"
)

# Italic formatting for pathogen names
italic_names <- c(
  "Acinetobacter spp."    = "<i>Acinetobacter</i> spp.",
  "Pseudomonas spp."      = "<i>Pseudomonas</i> spp.",
  "Enterococcus spp."     = "<i>Enterococcus</i> spp.",
  "Staphylococcus aureus" = "<i>Staphylococcus aureus</i>"
)

# Thresholds to skip under-powered subgroups
min_n_obs     <- 100     # minimum rows in subgroup
min_n_resist  <- 30      # minimum resistant cases in subgroup
min_n_events  <- 5       # minimum events for Cox PH

# Loop over each pathogen
for (j in seq_along(bug_labels)) {
  
  df_main    <- as.data.frame(df_all[[j]])
  df_used    <- split(df_main, df_main$infection_types)
  resist_var <- resist_vars[j]
  table_zph  <- list()
  
  # Loop over each infection-type subgroup
  for (i in seq_along(df_used)) {
    df  <- df_used[[i]]
    subgroup  <- names(df_used)[i]
    if (is.null(df) || nrow(df) == 0) {
      message("Skipped: empty subgroup '", subgroup, "'")
      next
    }
    
    # basic checks
    if (!(resist_var %in% names(df))) {
      message("Skipped: '", resist_var, "' not in subgroup '", subgroup, "'")
      next
    }
    
    n_obs    <- nrow(df)
    n_events <- sum(df$event, na.rm = TRUE)
    rv_chr   <- as.character(df[[resist_var]])
    n_resist <- sum(
      !is.na(rv_chr) &
        (rv_chr %in% c("R","r","1","resistant","Resistant","RESISTANT") |
           grepl("resistant$", rv_chr, ignore.case = TRUE)),
      na.rm = TRUE
    )
    
    # skip underpowered strata
    if (n_obs < min_n_obs) {
      message("Skipped: n_obs<", min_n_obs, " in subgroup '", subgroup, "' (n=", n_obs, ")")
      next
    }
    if (n_resist < min_n_resist) {
      message("Skipped: resistant<", min_n_resist, " in subgroup '", subgroup, "' (n_resist=", n_resist, ")")
      next
    }
    if (n_events < min_n_events) {
      message("Skipped: events<", min_n_events, " in subgroup '", subgroup, "' (events=", n_events, ")")
      next
    }
    
    # Build formula
    fml <- as.formula(paste(
      "Surv(time, event) ~ age_new + sex + country_region + country_income +",
      "hpd_admreason + comorbidities_CCI + severity_score_scale +",
      "icu_hd_ap + pathogen_combined_types +", resist_var
    ))

    
    # Fit Cox model with error catching
    fit_cox <- try(suppressWarnings(coxph(fml, data = df)), silent = TRUE)
    if (inherits(fit_cox, "try-error")) {
      message("Skipped: coxph() error in subgroup '", subgroup, "'")
      next
    }
    
    # Test PH assumption with error catching
    zph <- try(cox.zph(fit_cox), silent = TRUE)
    if (inherits(zph, "try-error")) {
      message("Skipped: cox.zph() error in subgroup '", subgroup, "'")
      next
    }
    
    # Format zph output
    df_zph       <- as.data.frame(zph$table)
    df_zph$p     <- ifelse(
      df_zph$p < 0.001, "<0.001",
      ifelse(df_zph$p < 0.01, sprintf("%.3f", df_zph$p),
             sprintf("%.2f", df_zph$p))
    )
    table_zph[[subgroup]] <- df_zph
    message("Processed ", bug_labels[j], " – ", subgroup,
            ": n=", n_obs, ", resist=", n_resist, ", events=", n_events)
  }
  
  # If all subgroups skipped, do not write a file
  if (length(table_zph) == 0) {
    message("Skip writing for: ", bug_labels[j], " (all infection-type subgroups skipped)")
    next
  }
  
  # Define row labels
  base_label <- sub(" \\(.*\\)", "", bug_labels[j])
  
  row_labels <- c(
    "Age", "Sex", "Region", "World Bank income status",
    "Primary admission reason", "Charlson comorbidity index",
    "Severity score of disease", "Admission to ICU/HD at enrollment",
    "Type of pathogens", base_label, "GLOBAL"
  )
  
  
  # Combine into one data.frame (only kept subgroups)
  combined <- data.frame(
    Variables   = row_labels,
    do.call(cbind, table_zph),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  
  subgroups  <- names(table_zph)
  stat_names <- c("chisq", "df", "p")
  colnames(combined) <- c(
    "Variables",
    unlist(lapply(subgroups, function(s) paste0(stat_names, "_", s)))
  )
  rownames(combined) <- NULL
  
  # Labels for columns
  labels <- rep(list(md("Chi-sq"), md("DF"), md("*p*")), times = length(subgroups))
  names(labels) <- unlist(lapply(subgroups, function(s) paste0(stat_names, "_", s)))
  
  # Build gt table
  gt_tbl <- combined %>%
    gt() %>%
    fmt_number(columns = matches("^chisq"), decimals = 3)
  
  # Add spanners for each infection-type (only kept ones)
  for (nm in subgroups) {
    gt_tbl <- gt_tbl %>%
      tab_spanner(
        label   = md(paste0("**<span style='color:black;'>", nm, "</span>**")),
        columns = matches(paste0("_", nm, "$"))
      )
  }
  
  # Continue styling
  gt_tbl <- gt_tbl %>%
    cols_label(.list = labels) %>%
    cols_width(
      1 ~ px(165),
      2:ncol(combined) ~ px(50)
    ) %>%
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_column_labels(everything())
    ) %>%
    tab_options(
      column_labels.border.top.color    = "black",
      column_labels.border.top.width    = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(2),
      table_body.hlines.color           = "white",
      table_body.hlines.width           = px(0),
      table.border.bottom.color         = "white",
      table.border.bottom.width         = px(0),
      data_row.padding                  = px(0)
    ) %>%
    tab_style(
      style     = cell_borders(sides = "bottom", color = "black", weight = px(2)),
      locations = cells_body(rows = nrow(combined))
    ) %>%
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_body(rows = nrow(combined))
    ) %>%
    cols_align(
      align   = "center",
      columns = 2:ncol(combined)
    ) %>%
    tab_style(
      style     = cell_text(font = "Times New Roman", size = px(10)),
      locations = list(
        cells_title(groups = "title"),
        cells_title(groups = "subtitle"),
        cells_column_labels(),
        cells_body()
      )
    ) %>%
    tab_style(
      style     = cell_text(font = "Times New Roman", size = px(10), weight = "bold"),
      locations = cells_column_spanners()
    ) %>%
    tab_source_note(
      source_note = "Abbreviations: Chi-sq = Chi-Squared statistic; DF = Degrees of freedom."
    ) %>%
    tab_style(
      style     = cell_text(font = "Times New Roman", size = px(10)),
      locations = cells_source_notes()
    ) %>%
    tab_style(
      style     = cell_text(align = "left", v_align = "middle"),
      locations = cells_column_labels(columns = 1)
    )
  
  # Italicize pathogen row if needed
  if (bug_labels[j] %in% names(italic_names)) {
    gt_tbl <- gt_tbl %>%
      text_transform(
        locations = cells_body(
          columns = Variables,
          rows    = Variables == bug_labels[j]
        ),
        fn = function(x) html(italic_names[[bug_labels[j]]])
      )
  }
  
  # Save
  gtsave(gt_tbl, filename = output_paths[j])
  message("Saved: ", output_paths[j], " | subgroups: ", paste(subgroups, collapse = ", "))
}
