# Clear 
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 rstpm2,
                 gt,
                 extrafont)
})

# Set working directory and load fonts
setwd("./")
loadfonts() 

# Read in the list of six datasets
df_all <- readRDS("data/att_first28.RData")

# Spline degrees of freedom
df_sp <- 2:5

# Six resistance variables
resist_vars <- c("aci_car",
                 "ent_thir",
                 "ent_car",
                 "pse_car",
                 "entc_van",
                 "sa_meth",
                 "mdr_gnb",
                 "mdr",
                 "amr")

# Paths for saving each HTML output
output_paths <- paste0("output/table/AIC_att28_",
                       sprintf("%02d", seq_along(resist_vars)),
                       "_death.html")

# Loop over each resistance variable
for (i in seq_along(resist_vars)) {
  # Extract the i-th data frame and variable name
  df  <- as.data.frame(df_all[[i]])
  var <- resist_vars[i]
  
  #Special handling for VRE
  if (var == "entc_van") {
    df$infection_types <- droplevels(factor(df$infection_types))
    df$infection_types <- relevel(df$infection_types, ref = "Hospital-acquired BSI")
  }
  
  # Compute AIC for each candidate df
  aic_rp <- vapply(df_sp, function(dfi) {
    # Dynamically build the model formula including the current var
    formula_text <- paste0(
      "Surv(time, event) ~ age_new + sex + country_region + country_income + hpd_admreason + comorbidities_CCI + severity_score_scale + icu_hd_ap + pathogen_combined_types + infection_types + ",
      var
    )
    mod <- stpm2(
      formula = as.formula(formula_text),
      data = df,
      df = dfi
    )
    AIC(mod)
  }, numeric(1))
  
  # Combine into a data frame and sort by AIC
  aic_df <- data.frame(
    DF  = df_sp,
    AIC = aic_rp
  ) %>% 
    arrange(AIC)
  
  # Create a styled gt table
  tbl <- aic_df %>%
    gt() %>%
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
      style = cell_borders(
        sides = "bottom",
        color = "black",
        weight = px(2)
      ),
      locations = cells_body(
        rows = nrow(aic_df)
      )
    ) %>%
    cols_align(
      align = "center",
      columns = 1:2
    ) %>%
    tab_style(
      style = cell_text(
        font = c("Times New Roman"),
        size = px(10)
      ),
      locations = list(
        cells_title(groups = "title"),
        cells_title(groups = "subtitle"),
        cells_column_labels(),
        cells_body()
      )
    ) %>%
  tab_source_note(
    source_note = html("Abbreviations: DF = Degrees of freedom,<br>AIC = Akaike Information Criterion.")
  ) %>%
    tab_style(
      style = cell_text(
        font = "Times New Roman",
        size = px(10)
      ),
      locations = cells_source_notes()
    ) 
  # Save to HTML
  gtsave(data = tbl, filename = output_paths[i])
}
