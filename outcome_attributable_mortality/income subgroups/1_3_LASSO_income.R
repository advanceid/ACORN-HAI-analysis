# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, car, gt, extrafont)
})

# Working directory and fonts
setwd("./")
loadfonts()

# Load data
df_all <- readRDS("data/att_first28_income.RData")

# Pathogen labels and resistance variable names
bug_labels <- c("Acinetobacter spp.", "Enterobacterales (3GCRE)", 
                "Enterobacterales (CRE)", "Pseudomonas spp.",
                "Enterococcus spp.", "Staphylococcus aureus", "MDR GNB")

resist_vars <- c("aci_car", "ent_thir", "ent_car",
                 "pse_car", "entc_van", "sa_meth", "mdr_gnb")

output_paths <- paste0("output/table/LASSO_att28_income_",
                       c("01_cra", "02_3gcre", "03_cre", 
                         "04_crp", "05_vre", "06_mrsa", "07_mdr_gnb"), ".html")

# Thresholds to skip strata
min_n_obs    <- 100
min_n_resist <- 30

# Safe extractor: return NA if vif vector too short
get_vif <- function(vif_vals, k) {
  if (length(vif_vals) < k) return(NA_real_)
  x <- suppressWarnings(as.numeric(vif_vals[k]))
  if (length(x) == 0) return(NA_real_)
  round(x, 3)
}

# Loop over pathogens
for (j in seq_along(bug_labels)) {
  
  df_main <- as.data.frame(df_all[[j]])
  
  # VRE: drop unused levels and set reference 
  if (bug_labels[j] == "Enterococcus spp.") {
    df_main$infection_types <- droplevels(factor(df_main$infection_types))
    df_main$infection_types <- stats::relevel(df_main$infection_types, ref = "Hospital-acquired BSI")
  }
  
  df_used <- split(df_main, df_main$country_income)  # group by income
  table_vif <- list()
  kept_groups <- character(0)
  
  if (bug_labels[j] == "Enterococcus spp.") {
    bold_rows <- c(1, 3, 6, 11, 13, 15, 18, 21, 24)
  } else {
    bold_rows <- c(1, 3, 6, 11, 13, 15, 18, 21, 25)
  }
  
  for (i in seq_along(df_used)) {
    df <- df_used[[i]]
    if (is.null(df) || nrow(df) == 0) next
    
    n_obs <- nrow(df)
    rv_chr <- as.character(df[[resist_vars[j]]])
    n_resist <- sum(
      !is.na(rv_chr) &
        (rv_chr %in% c("R","r","1","resistant","Resistant","RESISTANT") |
           grepl("resistant$", rv_chr, ignore.case = TRUE)),
      na.rm = TRUE
    )
    if (n_obs < min_n_obs)    next
    if (n_resist < min_n_resist) next
    
    x.factors <- model.matrix(
      reformulate(c("df$sex", "df$hpd_admreason", "df$icu_hd_ap",
                    "df$pathogen_combined_types", "df$infection_types",
                    paste0("df$", resist_vars[j])))
    )[, -1]
    
    x.factors <- as.matrix(data.frame(
      x.factors, df$age_new, df$comorbidities_CCI, df$severity_score_scale
    ))
    
    x.factors.df <- as.data.frame(x.factors)
    model <- lm(rnorm(nrow(x.factors.df)) ~ ., data = x.factors.df)
    vif_df <- as.data.frame(vif(model))
    vif_df$vif_values <- round(as.numeric(vif_df[, 1]), 3)
    
    # Susceptible/resistant labels for display only
    sus_label <- switch(j,
                        "Carbapenem-susceptible",                 
                        "Third-generation cephalosporin-susceptible",       
                        "Carbapenem-susceptible",                      
                        "Carbapenem-susceptible",
                        "Vancomycin-susceptible",
                        "Methicillin-susceptible")          
    
    res_label <- switch(j,
                        "Carbapenem-resistant",
                        "Third-generation cephalosporin-resistant",
                        "Carbapenem-resistant",
                        "Carbapenem-resistant",
                        "Vancomycin-resistant",
                        "Methicillin-resistant")
    
    # Build VIF table column (use get_vif to avoid out-of-bounds)
    # table rows
    if (bug_labels[j] == "Enterococcus spp.") {
      table_row <- list(
        c("Age", NA),
        c(NA, get_vif(vif_df$vif_values, 9)),
        c("Sex", NA),
        c("Female", "—"),
        c("Male", get_vif(vif_df$vif_values, 1)),
        c("Primary admission reason", NA),
        c("Infectious disease", "—"),
        c("Gastrointestinal disorder", get_vif(vif_df$vif_values, 2)),
        c("Pulmonary disease", get_vif(vif_df$vif_values, 3)),
        c("Others", get_vif(vif_df$vif_values, 4)),
        c("Charlson comorbidity index", NA),
        c(NA, get_vif(vif_df$vif_values, 10)),
        c("Severity score of disease", NA),
        c(NA, get_vif(vif_df$vif_values, 11)),
        c("Admission to ICU/HD at enrollment", NA),
        c("No", "—"),
        c("Yes", get_vif(vif_df$vif_values, 5)),
        c("Type of pathogens", NA),
        c("Monomicrobial", "—"),
        c("Polymicrobial", get_vif(vif_df$vif_values, 6)),
        c("Infection syndromes", NA),
        c("Hospital-acquired BSI", "—"),
        c("Healthcare-associated BSI", get_vif(vif_df$vif_values, 7)),
        c(sub(" \\(.*\\)", "", bug_labels[j]), NA),
        c(sus_label, "—"),
        c(res_label, get_vif(vif_df$vif_values, 8))
      )
    } else {
      table_row <- list(
        c("Age", NA),
        c(NA, get_vif(vif_df$vif_values, 10)),
        c("Sex", NA),
        c("Female", "—"),
        c("Male", get_vif(vif_df$vif_values, 1)),
        c("Primary admission reason", NA),
        c("Infectious disease", "—"),
        c("Gastrointestinal disorder", get_vif(vif_df$vif_values, 2)),
        c("Pulmonary disease", get_vif(vif_df$vif_values, 3)),
        c("Others", get_vif(vif_df$vif_values, 4)),
        c("Charlson comorbidity index", NA),
        c(NA, get_vif(vif_df$vif_values, 11)),
        c("Severity score of disease", NA),
        c(NA, get_vif(vif_df$vif_values, 12)),
        c("Admission to ICU/HD at enrollment", NA),
        c("No", "—"),
        c("Yes", get_vif(vif_df$vif_values, 5)),
        c("Type of pathogens", NA),
        c("Monomicrobial", "—"),
        c("Polymicrobial", get_vif(vif_df$vif_values, 6)),
        c("Infection syndromes", NA),
        c("VAP", "—"),
        c("Hospital-acquired BSI", get_vif(vif_df$vif_values, 7)),
        c("Healthcare-associated BSI", get_vif(vif_df$vif_values, 8)),
        c(sub(" \\(.*\\)", "", bug_labels[j]), NA),
        c(sus_label, "—"),
        c(res_label, get_vif(vif_df$vif_values, 9))
      )
    }
    
    col_df <- do.call(rbind, table_row) %>%
      as.data.frame() %>%
      mutate(across(1, ~ ifelse(.x == "NA", "", .x)))
    
    names(col_df) <- c("Variables", names(df_used)[i])  # name the kept income group
    table_vif[[length(table_vif) + 1]] <- col_df
    kept_groups <- c(kept_groups, names(df_used)[i])
  }
  
  # If all income groups were skipped, don't write a file
  if (length(table_vif) == 0) {
    message("Skip writing for: ", bug_labels[j], " (all income groups skipped)")
    next
  }
  
  # Combine only kept groups
  table_vif_all <- table_vif[[1]]
  if (length(table_vif) > 1) {
    for (k in 2:length(table_vif)) {
      table_vif_all <- cbind(table_vif_all, table_vif[[k]][, 2, drop = FALSE])
    }
  }
  colnames(table_vif_all) <- c("Variables", kept_groups)
  
  # Create gt table (auto-fit widths; avoids errors when only 1 group remains)
  width_varcol <- if (bug_labels[j] == "Enterobacterales (3GCRE)") px(185) else px(165)
  n_cols <- ncol(table_vif_all)  # >= 2
  
  lasso_table <- table_vif_all %>%
    gt() %>%
    sub_missing(columns = 1:n_cols, missing_text = "") %>%
    tab_style(style = list(cell_text(weight = "bold")),
              locations = cells_column_labels(everything()))
  
  if (n_cols == 2) {
    lasso_table <- lasso_table %>%
      cols_width(1 ~ width_varcol, 2 ~ px(110))
  } else if (n_cols == 3) {
    lasso_table <- lasso_table %>%
      cols_width(1 ~ width_varcol, 2 ~ px(110), 3 ~ px(110))
  } else {
    lasso_table <- lasso_table %>%
      cols_width(
        1 ~ width_varcol, 
        2 ~ px(70),
        tidyselect::all_of(3:n_cols) ~ px(120)
      )
  }
  
  lasso_table <- lasso_table %>%
    tab_options(
      column_labels.border.top.color = "black",
      column_labels.border.top.width = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(2),
      table_body.hlines.color = "white",
      table_body.hlines.width = px(0),
      table.border.bottom.color = "white",
      table.border.bottom.width = px(0),
      data_row.padding = px(0),
      source_notes.padding = px(0)
    ) %>%
    tab_style(style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
              locations = cells_body(rows = nrow(table_vif_all))) %>%
    cols_align(align = "center", columns = 2:n_cols) %>%
    tab_style(style = cell_text(align = "left", v_align = "middle"),
              locations = cells_column_labels(columns = 1)) %>%
    tab_spanner(label = md("**<span style = 'color:black;'>Variance inflation factor</span>**"),
                columns = 2:n_cols) %>%
    tab_style(style = cell_text(font = "Times New Roman", size = px(10), weight = "bold"),
              locations = cells_column_spanners(spanners = tidyselect::matches("Variance inflation factor"))) %>%
    tab_style(style = list(cell_text(weight = "bold")),
              locations = cells_body(columns = 1, rows = bold_rows)) %>%
    tab_style(style = cell_text(font = "Times New Roman", size = px(10)),
              locations = list(cells_title(groups = "title"),
                               cells_title(groups = "subtitle"),
                               cells_column_labels(),
                               cells_body()))
  
  # Scientific names: italic + bold
  if (bug_labels[j] == "Acinetobacter spp.") {
    lasso_table <- lasso_table %>%
      text_transform(
        locations = cells_body(columns = Variables, rows = Variables == "Acinetobacter spp."),
        fn = function(x) html("<b><i>Acinetobacter</i> spp.</b>")
      )
  }
  if (bug_labels[j] == "Pseudomonas spp.") {
    lasso_table <- lasso_table %>%
      text_transform(
        locations = cells_body(columns = Variables, rows = Variables == "Pseudomonas spp."),
        fn = function(x) html("<b><i>Pseudomonas</i> spp.</b>")
      )
  }
  if (bug_labels[j] == "Enterococcus spp.") {
    lasso_table <- lasso_table %>%
      text_transform(
        locations = cells_body(columns = Variables, rows = Variables == "Enterococcus spp."),
        fn = function(x) html("<b><i>Enterococcus</i> spp.</b>")
      )
  }
  if (bug_labels[j] == "Staphylococcus aureus") {
    lasso_table <- lasso_table %>%
      text_transform(
        locations = cells_body(columns = Variables, rows = Variables == "Staphylococcus aureus"),
        fn = function(x) html("<b><i>Staphylococcus aureus</i></b>")
      )
  }
  
  # Save table
  gtsave(lasso_table, filename = output_paths[j])
}

