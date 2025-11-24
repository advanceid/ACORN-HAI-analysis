# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, 
                 survival,
                 gt, 
                 gtsummary, 
                 extrafont)
})

# Set working directory and load fonts
setwd("./")
loadfonts()

# Load data
df_all <- readRDS("data/att_first28.RData")

# Define pathogen labels & output paths
bug_labels <- c("Acinetobacter spp.", "Enterobacterales (3GCRE)", 
                "Enterobacterales (CRE)", "Pseudomonas spp.",
                "Enterococcus spp.", "Staphylococcus aureus", 
                "MDR GNB", "MDR", "AMR")
resist_vars <- c("aci_car", "ent_thir", "ent_car", 
                 "pse_car", "entc_van", "sa_meth",
                 "mdr_gnb", "mdr", "amr")
output_paths <- paste0("output/table/cox_zph_att28_overall_",
                       c("01_cra", "02_3gcre", "03_cre", 
                         "04_crp", "05_vre", "06_mrsa",
                         "07_mdr_gnb", "08_mdr", "09_amr"), ".html")

italic_names <- c(
  "Acinetobacter spp." = "<i>Acinetobacter</i> spp.",
  "Pseudomonas spp." = "<i>Pseudomonas</i> spp.",
  "Enterococcus spp." = "<i>Enterococcus</i> spp.",
  "Staphylococcus aureus" = "<i>Staphylococcus aureus</i>"
)

zph_vars <- c("Age", "Sex", 
              "Region", "World Bank income status",
              "Primary admission reason",
              "Charlson comorbidity index",
              "Severity score of disease",
              "Admission to ICU/HD at enrollment",
              "Type of pathogens",
              "Infection syndromes",
              "BUG_PLACEHOLDER",
              "GLOBAL")

for (j in seq_along(bug_labels)) {
  df <- as.data.frame(df_all[[j]])
  resist_var <- resist_vars[j]
  
  #Special handling for VRE
  if (bug_labels[j] == "Enterococcus spp.") {
    df$infection_types <- droplevels(factor(df$infection_types))
    df$infection_types <- relevel(df$infection_types, ref = "Hospital-acquired BSI")
  }
  
  formula <- as.formula(
    paste("Surv(time, event) ~ age_new + sex + country_region + country_income + hpd_admreason + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + infection_types +", resist_var)
  )
  
  model_cox <- suppressWarnings(coxph(formula, data = df))
  cox_zph_result <- cox.zph(model_cox)
  cox_zph_df <- as.data.frame(cox_zph_result$table)
  
  # Format p-values
  cox_zph_df$p <- ifelse(cox_zph_df$p < 0.001, "<0.001",
                         ifelse(cox_zph_df$p < 0.01, sprintf("%.3f", cox_zph_df$p), 
                                sprintf("%.2f", cox_zph_df$p)))
  
  table_zph <- data.frame(
    Variables = sub("BUG_PLACEHOLDER", sub(" \\(.*\\)", "", bug_labels[j]), zph_vars),
    Chi_sq = round(cox_zph_df$chisq, 3),
    DF = cox_zph_df$df,
    p = cox_zph_df$p,
    check.names = FALSE
  )
  
  # Build gt table
  table <- table_zph %>%
    gt() %>%
    cols_label(
      Chi_sq = md("Chi-sq"),
      DF = md("DF"),
      p = md("*p*")
    ) %>%
    cols_width(
      Variables ~ px(165),
      everything() ~ px(50)
    ) %>%
    fmt_number(columns = Chi_sq, decimals = 3) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels(everything())
    ) %>%
    tab_style(
      style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
      locations = cells_body(rows = nrow(table_zph))
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
    cols_align(
      align = "center",
      columns = c(Chi_sq, DF, p)
    ) %>%
    tab_style(
      style = cell_text(align = "left", v_align = "middle"),
      locations = cells_column_labels(columns = Variables)
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
    tab_source_note(
      source_note = "Abbreviations: Chi-sq = Chi-Squared Statistic, DF = Degrees of freedom."
    ) %>%
    tab_style(
      style = cell_text(font = "Times New Roman", size = px(10)),
      locations = cells_source_notes()
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(rows = Variables == "GLOBAL")
    )
  
  # Italicize if needed
  if (bug_labels[j] %in% names(italic_names)) {
    table <- table %>%
      text_transform(
        locations = cells_body(columns = Variables, rows = Variables == bug_labels[j]),
        fn = function(x) html(italic_names[[bug_labels[j]]])
      )
  }
  
  gtsave(data = table, filename = output_paths[j])
}