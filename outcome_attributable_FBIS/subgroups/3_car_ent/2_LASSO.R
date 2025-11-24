# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 ordinalNet,
                 car,
                 gt,
                 extrafont)
})

# Working directory
wd <- "./"
setwd(wd)

# Load dataset
df_all <- readRDS("data/att_fbis.RData")

# Load fonts
loadfonts()

# Thresholds and safe VIF extractor
min_n_obs    <- 100   # minimum number of observations
min_n_resist <- 30    # minimum number of resistant cases

get_vif <- function(vif_vals, k) {
  if (length(vif_vals) < k) return(NA_real_)
  x <- suppressWarnings(as.numeric(vif_vals[k]))
  if (length(x) == 0) return(NA_real_)
  round(x, 3)
}

# Enterobacterales (CRE) — using ent_car as resistance var
df_all <- as.data.frame(df_all[[3]])
df_used <- split(df_all, df_all$infection_types)

table_vif   <- list()
kept_groups <- character(0)

system.time({
  for (i in seq_along(df_used)) {
    df <- df_used[[i]]
    if (is.null(df) || nrow(df) == 0) next
    
    # Skip subgroups with too few observations/resistant cases
    n_obs <- nrow(df)
    rv_chr <- as.character(df$ent_car)
    n_resist <- sum(
      !is.na(rv_chr) &
        (rv_chr %in% c("R","r","1","resistant","Resistant","RESISTANT") |
           grepl("resistant$", rv_chr, ignore.case = TRUE)),
      na.rm = TRUE
    )
    if (n_obs < min_n_obs)       next
    if (n_resist < min_n_resist) next
    
    # Design matrix
    x.factors <- model.matrix(
      ~ df$sex + df$country_region +
        df$country_income + df$hpd_admreason +
        df$icu_hd_ap + df$pathogen_combined_types +
        df$ent_car
    )[,-1]
    
    # Add numeric covariates
    x.factors <- as.matrix(data.frame(
      x.factors,
      df$age_new,
      df$comorbidities_CCI,
      df$severity_score_scale
    ))
    
    # Dummy outcome just for VIF calculation
    x.factors.df <- data.frame(x.factors)
    dummy_y <- rnorm(nrow(x.factors.df))
    model <- lm(dummy_y ~ ., data = x.factors.df)
    
    # Extract VIF values
    vif_values <- car::vif(model)
    vif_df  <- as.data.frame(vif_values)
    vif_vals <- round(as.numeric(vif_df[,1]), 3)
    
    # For table
    table_row <- list(
      c("Age", NA),
      c(NA, get_vif(vif_vals, 12)),
      c("Sex", NA),
      c("Female", "—"),
      c("Male", get_vif(vif_vals, 1)),
      c("Region", NA),
      c("Eastern Mediterranean Region", "—"),
      c("South-East Asian Region", get_vif(vif_vals, 2)),
      c("Western Pacific Region", get_vif(vif_vals, 3)),
      c("World Bank income status", NA),
      c("High income", "—"),
      c("Upper middle income", get_vif(vif_vals, 4)),
      c("Lower middle income", get_vif(vif_vals, 5)),
      c("Primary admission reason", NA),
      c("Infectious disease", "—"),
      c("Gastrointestinal disorder", get_vif(vif_vals, 6)),
      c("Pulmonary disease", get_vif(vif_vals, 7)),
      c("Others", get_vif(vif_vals, 8)),
      c("Charlson comorbidity index", NA),
      c(NA, get_vif(vif_vals, 13)),
      c("Severity score of disease", NA),
      c(NA, get_vif(vif_vals, 14)),
      c("Admission to ICU/HD at enrollment", NA),
      c("No", "—"),
      c("Yes", get_vif(vif_vals, 9)),
      c("Type of pathogens", NA),
      c("Monomicrobial", "—"),
      c("Polymicrobial", get_vif(vif_vals, 10)),
      c("Enterobacterales", NA),
      c("Carbapenem-susceptible", "—"),
      c("Carbapenem-resistant", get_vif(vif_vals, 11))
    )
    
    col_df <- do.call(rbind, table_row) |>
      as.data.frame() |>
      dplyr::mutate(across(1, ~ ifelse(.x == "NA", "", .x)))
    names(col_df) <- c("Variables", names(df_used)[i])
    
    table_vif[[length(table_vif) + 1]] <- col_df
    kept_groups <- c(kept_groups, names(df_used)[i])
  }
})

# Stop if all subgroups skipped
if (length(table_vif) == 0) {
  message("All subgroups skipped (below thresholds). No file written.")
  quit(save = "no")
}

# Merge only retained subgroups
table_vif_all <- table_vif[[1]]
if (length(table_vif) > 1) {
  for (k in 2:length(table_vif)) {
    table_vif_all <- cbind(table_vif_all, table_vif[[k]][, 2, drop = FALSE])
  }
}
colnames(table_vif_all) <- c("Variables", kept_groups)

# gt table
n_cols <- ncol(table_vif_all)

lasso_table <- table_vif_all %>%
  gt() %>%
  sub_missing(columns = 1:n_cols, missing_text = "") %>%
  tab_style(style = list(cell_text(weight = "bold")),
            locations = cells_column_labels(everything()))

if (n_cols == 2) {
  lasso_table <- lasso_table %>%
    cols_width(1 ~ px(165), 2 ~ px(110))
} else if (n_cols == 3) {
  lasso_table <- lasso_table %>%
    cols_width(1 ~ px(165), 2 ~ px(80), 3 ~ px(120))
} else {
  lasso_table <- lasso_table %>%
    cols_width(1 ~ px(165), 2 ~ px(60), tidyselect::all_of(3:n_cols) ~ px(125))
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
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
    locations = cells_body(rows = nrow(table_vif_all))
  ) %>%
  cols_align(align = "center", columns = 2:n_cols) %>%
  tab_style(style = cell_text(align = "left", v_align = "middle"),
            locations = cells_column_labels(columns = 1)) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Variance inflation factor</span>**"),
    columns = 2:n_cols
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman", size = px(10), weight = "bold"),
    locations = cells_column_spanners(spanners = tidyselect::matches("Variance inflation factor"))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = 1,
      rows = c(1, 3, 6, 10, 14, 19, 21, 23, 26, 29)
    )
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman", size = px(10)),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  )

# Preview
lasso_table

# Save
gtsave(lasso_table, filename = "output/table/LASSO_attfbis.html")
