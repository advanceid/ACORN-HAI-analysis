# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 glmnet,
                 car,
                 gt,
                 extrafont)
})
# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_fbis.RData")

# Load fonts
loadfonts()

# Enterobacterales, third-generation cephalosporin-resistant
df <- as.data.frame(df_all[[2]])

# LASSO
# For variable selection when there is multicollinearity among the predictors
# Matrix for factor
x.factors <- model.matrix(~ df$sex + df$country_region +
                            df$country_income + 
                            df$hpd_admreason +
                            df$icu_hd_ap + df$pathogen_combined_types+
                            df$ent_thir)[,-1]

# For num
x.factors <- as.matrix(data.frame(x.factors, df$age_new, 
                                  df$comorbidities_CCI, 
                                  df$severity_score_scale))
# For factor
x.factors.df <- data.frame(x.factors)
dummy_y <- rnorm(nrow(x.factors.df))

model <- lm(dummy_y ~ ., data = x.factors.df)

# Variance Inflation Factor (VIF)
vif_values <- vif(model)
vif_df <- as.data.frame(vif_values)
vif_df$vif_values <- round(vif_df$vif_values, 3)

# Creat table
table_row <- list(
  c("Age", NA),
  c(NA, vif_df$vif_values[12]),
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
  c(NA, vif_df$vif_values[13]),
  c("Severity score of disease", NA),
  c(NA, vif_df$vif_values[14]),
  c("Admission to ICU/HD at enrollment", NA),
  c("No", "—"),
  c("Yes", vif_df$vif_values[9]),
  c("Type of pathogens", NA),
  c("Monomicrobial", "—"),
  c("Polymicrobial", vif_df$vif_values[10]),
  c("Enterobacterales", NA),
  c("Third-generation cephalosporin-susceptible", "—"),
  c("Third-generation cephalosporin-resistant", vif_df$vif_values[11])
)

# Bind the rows into a matrix
table_vif <- do.call(rbind, table_row)

#
table_vif <- table_vif %>%
  as.data.frame() %>%
  mutate(across(1, ~ ifelse(.x == "NA", "", .x)))

colnames(table_vif) <- c("Variables", "Variance inflation factor")

#
lasso_table <- table_vif %>%
  gt() %>%
  sub_missing(
    columns = 1:2,
    missing_text = ""
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
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(
      rows = nrow(table_vif)
    )
  ) %>%
  cols_align(
    align = "center",
    columns = 2
  )  %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1)
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = 1,
      rows = c(1, 3, 6, 10, 14, 19, 21, 23, 26, 29) 
    )
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
  ) 

print(lasso_table)

# Save
gtsave(lasso_table, filename = "output/table/LASSO_attfbis_overall_2.html")
