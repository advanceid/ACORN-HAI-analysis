# Clear
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 tidyr,
                 magrittr,
                 gt,
                 gtsummary,
                 extrafont,
                 openxlsx,
                 purrr)
})

# Set working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df_res1 <- readRDS("data/excess_los_overall.RData")
df_res2 <- readRDS("data/excess_los_by_infection_type.RData")
df_res3 <- readRDS("data/excess_los_by_income.RData")

# Pivot each to wide format
res1_wide <- df_res1 %>%
  filter(Cutoff %in% c("t14","t21","t28","t60","t90")) %>%
  select(group_by, Group, Pattern, Cutoff, Estimate_CI) %>%
  pivot_wider(
    id_cols     = c(group_by, Group, Pattern),
    names_from  = Cutoff,
    values_from = Estimate_CI
  )

res2_wide <- df_res2 %>%
  filter(Cutoff %in% c("t14","t21","t28","t60","t90")) %>%
  select(group_by, Group, Pattern, Cutoff, Estimate_CI) %>%
  pivot_wider(
    id_cols     = c(group_by, Group, Pattern),
    names_from  = Cutoff,
    values_from = Estimate_CI
  )

res3_wide <- df_res3 %>%
  filter(Cutoff %in% c("t14","t21","t28","t60","t90")) %>%
  select(group_by, Group, Pattern, Cutoff, Estimate_CI) %>%
  pivot_wider(
    id_cols     = c(group_by, Group, Pattern),
    names_from  = Cutoff,
    values_from = Estimate_CI
  )

# Function for combine
section_block <- function(title, rows1 = NULL, rows2 = NULL, rows3 = NULL) {
  # 1) header row
  header <- tibble(
    `Priority pathogens` = title,
    Group = NA_character_,
    t14 = NA_character_,
    t21 = NA_character_,
    t28 = NA_character_,
    t60 = NA_character_,
    t90 = NA_character_
  )
  
  # 2) body rows
  body1 <- res1_wide %>%
    slice(rows1) %>%
    transmute(
      `Priority pathogens` = group_by,
      Group,
      t14, t21, t28, t60, t90
    )
  
  body2 <- res2_wide %>%
    slice(rows2) %>%
    transmute(
      `Priority pathogens` = group_by,
      Group,
      t14, t21, t28, t60, t90
    )
  
  body3 <- res3_wide %>%
    slice(rows3) %>%
    transmute(
      `Priority pathogens` = group_by,
      Group,
      t14, t21, t28, t60, t90
    )
  
  # 3) stack them
  bind_rows(header, body1, body2, body3)
}

# List of blocks
blocks <- list(
  section_block("Carbapenem-resistant Acinetobacter spp.",  1, 1:2, 1:2),
  section_block("Third-generation cephalosporin-resistant Enterobacterales", 2, 4, 3:5),
  section_block("Carbapenem-resistant Enterobacterales", 3, 6:7, 6:8)
)

# Combine
result <- bind_rows(blocks) %>%
  replace(is.na(.), "") %>%
  mutate(across(where(is.character),
                ~ gsub("(?<![0-9])-0(?![0-9])", "0", ., perl = TRUE)))

result_90 <- result[,c(1,2,7)]


colnames(result_90) <- c("Critical pathogens", " ", "Excess LOS (95% CI), days")

#
table <- result_90 %>%
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
    data_row.padding = px(0),
    footnotes.padding = px(0),
    source_notes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(
      rows = nrow(result_90)
    )
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      rows = c(1, 7, 13)
    )
  )%>%
  cols_align(
    align = "center",
    columns = 2:ncol(result_90)
  ) %>%
  cols_align(
    align = "left",
    columns = 1
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
  tab_style(
    style = cell_text(
      font = c("Times New Roman"),
      size = px(10)
    ),
    locations = cells_footnotes()
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10)
    ),
    locations = cells_source_notes()
  ) %>%
  cols_width(
    1 ~ px(270),
    2 ~ px(160),
    3 ~ px(200)
  ) %>%
  tab_source_note(
    source_note = c("Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, LOS = Length Of Stay, CI = Confidence Interval.")
  ) %>%
  tab_style(
    style = cell_text(align = "center", v_align = "middle"),
    locations = cells_column_labels(columns = 1:ncol(result_90))
  ) %>%
  text_transform(
    locations = cells_body(
      columns = `Critical pathogens`,
      rows = `Critical pathogens` == "Carbapenem-resistant Acinetobacter spp."
    ),
    fn = function(x) html("<b>Carbapenem-resistant <i>Acinetobacter</i> spp.</b>")
  ) %>%
  tab_footnote(
    footnote = "Subgroups with total n < 100 or resistant cases < 30 were excluded.",
    locations = cells_body(
      columns = "Critical pathogens",
      rows    = `Critical pathogens` %in% 
        c("Infection syndromes", "World Bank income status")
    )
  ) %>%
  tab_options(footnotes.marks = c("‡"))

print(table)

# Save
gtsave(data = table, filename = "excess_LOS_90.html")
write.xlsx(result, "excess_LOS_90.xlsx", rowNames = FALSE)