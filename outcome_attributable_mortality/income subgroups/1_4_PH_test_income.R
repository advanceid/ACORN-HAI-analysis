# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, survival, gt, gtsummary, extrafont)
})

# Set working directory and load fonts
setwd("./")
loadfonts()

# Load data
df_all <- readRDS("data/att_first28_income.RData")

# Pathogen labels & output paths
bug_labels <- c("Acinetobacter spp.", "Enterobacterales (3GCRE)", 
                "Enterobacterales (CRE)", "Pseudomonas spp.",
                "Enterococcus spp.", "Staphylococcus aureus", "MDR GNB")
resist_vars <- c("aci_car", "ent_thir", "ent_car", 
                 "pse_car", "entc_van", "sa_meth", "mdr_gnb")
output_paths <- paste0("output/table/cox_zph_att28_income_",
                       c("01_cra", "02_3gcre", "03_cre", 
                         "04_crp", "05_vre", "06_mrsa", "07_mdr"), ".html")

# Italic formatting for scientific names
italic_names <- c(
  "Acinetobacter spp." = "<i>Acinetobacter</i> spp.",
  "Pseudomonas spp." = "<i>Pseudomonas</i> spp.",
  "Enterococcus spp." = "<i>Enterococcus</i> spp.",
  "Staphylococcus aureus" = "<i>Staphylococcus aureus</i>"
)

# Variable labels
zph_vars <- c("Age", "Sex", 
              "Primary admission reason",
              "Charlson comorbidity index",
              "Severity score of disease",
              "Admission to ICU/HD at enrollment",
              "Type of pathogens",
              "Infection syndromes",
              "BUG_PLACEHOLDER",
              "GLOBAL")

# ---- skip thresholds ----
min_n_obs    <- 100
min_n_resist <- 30

# Loop through pathogens
for (j in seq_along(bug_labels)) {
  
  df_main <- as.data.frame(df_all[[j]])
  
  # VRE: drop unused levels and set reference 
  if (bug_labels[j] == "Enterococcus spp.") {
    df_main$infection_types <- droplevels(factor(df_main$infection_types))
    df_main$infection_types <- stats::relevel(df_main$infection_types, ref = "Hospital-acquired BSI")
  }
  
  df_used <- split(df_main, df_main$country_income)
  table_zph <- list()
  kept_groups <- character(0)
  resist_var <- resist_vars[j]
  
  for (i in seq_along(df_used)) {
    df <- df_used[[i]]
    if (is.null(df) || nrow(df) == 0) next
    if (!(resist_var %in% names(df))) next
    
    # Skip rule
    n_obs <- nrow(df)
    rv_chr <- as.character(df[[resist_var]])
    n_resist <- sum(
      !is.na(rv_chr) &
        (rv_chr %in% c("R","r","1","resistant","Resistant","RESISTANT") |
           grepl("resistant$", rv_chr, ignore.case = TRUE)),
      na.rm = TRUE
    )
    if (n_obs < min_n_obs)    next
    if (n_resist < min_n_resist) next
    
    # Cox model
    formula <- as.formula(
      paste("Surv(time, event) ~ age_new + sex + hpd_admreason + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + infection_types +",
            resist_var)
    )
    model_cox <- suppressWarnings(coxph(formula, data = df))
    cox_zph_result <- cox.zph(model_cox)
    cox_zph_df <- as.data.frame(cox_zph_result$table)
    
    # format p-values
    cox_zph_df$p <- ifelse(cox_zph_df$p < 0.001, "<0.001",
                           ifelse(cox_zph_df$p < 0.01, sprintf("%.3f", cox_zph_df$p), 
                                  sprintf("%.2f", cox_zph_df$p)))
    table_zph[[length(table_zph) + 1]] <- cox_zph_df
    kept_groups <- c(kept_groups, names(df_used)[i])
  }
  
  # If all groups skipped, no file output
  if (length(table_zph) == 0) {
    message("Skip writing for: ", bug_labels[j], " (all income groups skipped)")
    next
  }
  
  # Combine results
  table_zph_all <- cbind(
    Variables = sub("BUG_PLACEHOLDER", sub(" \\(.*\\)", "", bug_labels[j]), zph_vars),
    do.call(cbind, table_zph)
  )
  
  colnames(table_zph_all) <- c("Variables",
                               paste0(rep(c("chisq", "df", "p"), times = length(kept_groups)),
                                      "_", rep(kept_groups, each = 3)))
  rownames(table_zph_all) <- NULL
  
  # Build gt table
  table <- table_zph_all %>%
    gt() %>%
    fmt_number(columns = matches("^chisq"), decimals = 3) %>%
    # Dynamically build spanners for kept_groups
    {
      tab_obj <- .
      for (grp in kept_groups) {
        label <- paste0("**<span style = 'color:black;'>", grp, "</span>**")
        tab_obj <- tab_obj %>%
          tab_spanner(
            label = md(label),
            columns = matches(paste0("_", grp, "$")),
            id = grp
          )
      }
      tab_obj
    } %>%
    cols_label(.list = setNames(rep(list(md("Chi-sq"), md("DF"), md("*p*")), length(kept_groups)),
                                paste0(rep(c("chisq", "df", "p"), times = length(kept_groups)),
                                       "_", rep(kept_groups, each = 3)))) %>%
    cols_width(
      1 ~ px(165),
      2:ncol(table_zph_all) ~ px(50)
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_column_labels(everything())
    ) %>%
    tab_options(
      column_labels.border.top.color = "black",
      column_labels.border.top.width = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(2),
      table_body.hlines.color = "white",
      table_body.hlines.width = px(0),
      table.border.bottom.color = "white",
      table.border.bottom.width = px(0),
      data_row.padding = px(0)
    ) %>%
    tab_style(
      style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
      locations = cells_body(rows = nrow(table_zph_all))
    ) %>%
    cols_align(
      align = "center",
      columns = 2:ncol(table_zph_all)
    ) %>%
    tab_style(
      style = cell_text(font = "Times New Roman", size = px(10)),
      locations = list(
        cells_title(groups = "title"),
        cells_title(groups = "subtitle"),
        cells_column_labels(),
        cells_body()
      )
    ) %>%
    tab_style(
      style = cell_text(font = "Times New Roman", size = px(10), weight = "bold"),
      locations = cells_column_spanners(spanners = kept_groups)
    ) %>%
    tab_source_note(
      source_note = "Abbreviations: Chi-sq = Chi-Squared Statistic, DF = Degrees of freedom."
    ) %>%
    tab_style(
      style = cell_text(font = "Times New Roman", size = px(10)),
      locations = cells_source_notes()
    ) %>%
    tab_style(
      style = cell_text(align = "left", v_align = "middle"),
      locations = cells_column_labels(columns = 1)
    )
  
  # Italicize scientific names
  if (bug_labels[j] %in% names(italic_names)) {
    table <- table %>%
      text_transform(
        locations = cells_body(columns = Variables, rows = Variables == bug_labels[j]),
        fn = function(x) html(italic_names[[bug_labels[j]]])
      )
  }
  
  # Save output
  gtsave(data = table, filename = output_paths[j])
  message("Saved: ", output_paths[j], " | groups: ", paste(kept_groups, collapse = ", "))
}
