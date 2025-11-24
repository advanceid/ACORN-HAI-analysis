# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, car, gt, extrafont)
})

# Set working directory and load fonts
setwd("./")
loadfonts()

# Load data
df_all <- readRDS("data/att_first28.RData")

# Define pathogen labels and resistance variable names
bug_labels <- c("Acinetobacter spp.", "Enterobacterales (3GCRE)", 
                "Enterobacterales (CRE)", "Pseudomonas spp.",
                "Enterococcus spp.", "Staphylococcus aureus", 
                "MDR GNB", "MDR", "AMR")
resist_vars <- c("aci_car", "ent_thir", "ent_car",
                 "pse_car", "entc_van", "sa_meth", 
                 "mdr_gnb", "mdr", "amr")
output_paths <- paste0("output/table/LASSO_att28_overall_",
                       c("01_cra", "02_3gcre", "03_cre", 
                         "04_crp", "05_vre", "06_mrsa",
                         "07_mdr_gnb", "08_mdr", "09_amr"), ".html")


for (j in seq_along(bug_labels)) {
  df <- as.data.frame(df_all[[j]])
  
  # VRE: drop unused levels and set reference 
  if (bug_labels[j] == "Enterococcus spp.") {
    df$infection_types <- droplevels(factor(df$infection_types))
    df$infection_types <- stats::relevel(df$infection_types, ref = "Hospital-acquired BSI")
  }
  
  if (bug_labels[j] == "Enterococcus spp.") {
    bold_rows <- c(1, 3, 6, 10, 14, 19, 21, 23, 26, 29, 32)
  } else {
    bold_rows <- c(1, 3, 6, 10, 14, 19, 21, 23, 26, 29, 33)
  }
  
  
  # model matrix
  x.factors <- model.matrix(
    reformulate(c("df$sex", "df$country_region","df$country_income", 
                  "df$hpd_admreason", "df$icu_hd_ap",
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
  
  sus_label <- switch(j,
                      "Carbapenem-susceptible",                 
                      "Third-generation cephalosporin-susceptible",       
                      "Carbapenem-susceptible",                      
                      "Carbapenem-susceptible",
                      "Vancomycin-susceptible",
                      "Methicillin-susceptible",
                      "Non-MDR GNB",
                      "Non-MDR",
                      "Non-AMR")          
  
  res_label <- switch(j,
                      "Carbapenem-resistant",
                      "Third-generation cephalosporin-resistant",
                      "Carbapenem-resistant",
                      "Carbapenem-resistant",
                      "Vancomycin-resistant",
                      "Methicillin-resistant",
                      "MDR GNB",
                      "MDR",
                      "AMR")
  
  # table rows
  if (bug_labels[j] == "Enterococcus spp.") {
    table_row <- list(
      c("Age", NA),
      c(NA, vif_df$vif_values[13]),
      c("Sex", NA),
      c("Female", "—"),
      c("Male", vif_df$vif_values[1]),
      c("Region", NA),
      c("Eastern Mediterranean Region", "—"),
      c("South-East Asian Region", vif_df$vif_values[2]),
      c("Western Pacific Region", vif_df$vif_values[3]),
      c("World Bank income status", NA),
      c("High income", "—"),
      c("Upper middle income", vif_df$vif_values[4]),
      c("Lower middle income", vif_df$vif_values[5]),
      c("Primary admission reason", NA),
      c("Infectious disease", "—"),
      c("Gastrointestinal disorder", vif_df$vif_values[6]),
      c("Pulmonary disease", vif_df$vif_values[7]),
      c("Others", vif_df$vif_values[8]),
      c("Charlson comorbidity index", NA),
      c(NA, vif_df$vif_values[14]),
      c("Severity score of disease", NA),
      c(NA, vif_df$vif_values[15]),
      c("Admission to ICU/HD at enrollment", NA),
      c("No", "—"),
      c("Yes", vif_df$vif_values[9]),
      c("Type of pathogens", NA),
      c("Monomicrobial", "—"),
      c("Polymicrobial", vif_df$vif_values[10]),
      c("Infection syndromes", NA),
      c("Hospital-acquired BSI", "—"),    
      c("Healthcare-associated BSI", vif_df$vif_values[11]),
      c(sub(" \\(.*\\)", "", bug_labels[j]), NA),
      c(sus_label, "—"),
      c(res_label, vif_df$vif_values[12])
    )
  } else {
    table_row <- list(
      c("Age", NA),
      c(NA, vif_df$vif_values[14]),
      c("Sex", NA),
      c("Female", "—"),
      c("Male", vif_df$vif_values[1]),
      c("Region", NA),
      c("Eastern Mediterranean Region", "—"),
      c("South-East Asian Region", vif_df$vif_values[2]),
      c("Western Pacific Region", vif_df$vif_values[3]),
      c("World Bank income status", NA),
      c("High income", "—"),
      c("Upper middle income", vif_df$vif_values[4]),
      c("Lower middle income", vif_df$vif_values[5]),
      c("Primary admission reason", NA),
      c("Infectious disease", "—"),
      c("Gastrointestinal disorder", vif_df$vif_values[6]),
      c("Pulmonary disease", vif_df$vif_values[7]),
      c("Others", vif_df$vif_values[8]),
      c("Charlson comorbidity index", NA),
      c(NA, vif_df$vif_values[15]),
      c("Severity score of disease", NA),
      c(NA, vif_df$vif_values[16]),
      c("Admission to ICU/HD at enrollment", NA),
      c("No", "—"),
      c("Yes", vif_df$vif_values[9]),
      c("Type of pathogens", NA),
      c("Monomicrobial", "—"),
      c("Polymicrobial", vif_df$vif_values[10]),
      c("Infection syndromes", NA),
      c("VAP", "—"),
      c("Hospital-acquired BSI", vif_df$vif_values[11]),
      c("Healthcare-associated BSI", vif_df$vif_values[12]),
      c(sub(" \\(.*\\)", "", bug_labels[j]), NA),
      c(sus_label, "—"),
      c(res_label, vif_df$vif_values[13])
    )
  }

  table_vif <- do.call(rbind, table_row) %>%
    as.data.frame() %>%
    mutate(across(1, ~ ifelse(.x == "NA", "", .x)))
  
  colnames(table_vif) <- c("Variables", "Variance inflation factor")
  
  lasso_table <- table_vif %>%
    gt() %>%
    sub_missing(missing_text = "") %>%
    cols_width(1 ~ px(185), 2 ~ px(120)) %>%
    cols_align("center", columns = 2) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(everything())) %>%
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
              locations = cells_body(rows = nrow(table_vif))) %>%
    tab_style(style = list(cell_text(weight = "bold")),
              locations = cells_body(columns = 1, rows = bold_rows)) %>%
    tab_style(style = cell_text(font = "Times New Roman", size = px(10)),
              locations = list(cells_title(groups = "title"),
                               cells_title(groups = "subtitle"),
                               cells_column_labels(),
                               cells_body()))
  
  # Italic + bold scientific names
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
  
  gtsave(lasso_table, filename = output_paths[j])
}