# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 gt,
                 purrr,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Combine
res1 <- readRDS("data/attfbis_table_multi_1_sub.RData")
res2 <- readRDS("data/attfbis_table_multi_2_sub.RData")

#
res_vap <- cbind.data.frame(res1[[1]][,c(1,5,6,10)], 
                            res2[[1]][,c(5,6,10)])

res_bsi_1 <- cbind.data.frame(res1[[2]][,c(5,6,10)], 
                              res2[[2]][,c(5,6,10)])

res <- cbind.data.frame(res_vap, res_bsi_1)

# Set col names
col_names <- as.character(unlist(res[1, ]))
res <- res[-1, ]
rownames(res) <- NULL

colnames(res) <- c("Characteristics", 
                   "Unstandardized coefficient (95%CI)_vap_1", 
                   "P value_vap_1",
                   "Standardized coefficient (95%CI)_vap_1", 
                   "Unstandardized coefficient (95%CI)_vap_2", 
                   "P value_vap_2",
                   "Standardized coefficient (95%CI)_vap_2",
                   "Unstandardized coefficient (95%CI)_bsi1_1", 
                   "P value_bsi1_1",
                   "Standardized coefficient (95%CI)_bsi1_1", 
                   "Unstandardized coefficient (95%CI)_bsi1_2", 
                   "P value_bsi1_2",
                   "Standardized coefficient (95%CI)_bsi1_2")
# 
process_missing_values <- function(...) {
  row <- list(...)
  for (i in 2:ncol(res)) {
    if (is.na(row[[i]])) {
      if (i > 1 && (is.na(row[[i-1]]) || row[[i-1]] == "")) {
        row[[i]] <- ""
      } else if (i > 1 && row[[i-1]] == "Ref") {
        row[[i]] <- "—"
      }
    }
  }
  return(row)
}

#
res <- res %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x) | .x == "NA", "", .x))) %>%
  pmap_dfr(process_missing_values) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "", .x))) 

# Create your gt table
table <- res %>%
  gt() %>%
  sub_missing(
    columns = 1:ncol(res),
    missing_text = ""
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels(everything())
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>VAP</span></sup>**"),
    columns = 2:7,
    level = 2,
    id = "vap"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup> </sup>**"),
    columns = 2:4,
    level = 1,
    id = "vap_sem_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup>[1]</sup>**"),
    columns = 5:7,
    level = 1,
    id = "vap_sem_2"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Hospital-acquired BSI</span></sup>**"),
    columns = 8:13,
    level = 2,
    id = "bsi1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup> </sup>**"),
    columns = 8:10,
    level = 1,
    id = "bsi1_sem_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup>[1]</sup>**"),
    columns = 11:13,
    level = 1,
    id = "bsi1_sem_2"
  ) %>%
  cols_label(
    `Unstandardized coefficient (95%CI)_vap_1` = md("Unstandardized coefficient (95%CI))"),
    `P value_vap_1` = md("*p*"),
    `Standardized coefficient (95%CI)_vap_1` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_vap_2` = md("Unstandardized coefficient (95%CI)"),
    `P value_vap_2` = md("*p*"),
    `Standardized coefficient (95%CI)_vap_2` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_bsi1_1` = md("Unstandardized coefficient (95%CI))"),
    `P value_bsi1_1` = md("*p*"),
    `Standardized coefficient (95%CI)_bsi1_1` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_bsi1_2` = md("Unstandardized coefficient (95%CI)"),
    `P value_bsi1_2` = md("*p*"),
    `Standardized coefficient (95%CI)_bsi1_2` = md("Standardized coefficient (95%CI)")
  ) %>%
  cols_align(
    align = "center",
    columns = 2:ncol(res)
  ) %>%
  cols_align(
    align = "left",
    columns = 1
  ) %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1)
  ) %>%
  tab_style(
    style = cell_text(align = "center", v_align = "middle"),
    locations = cells_column_labels(columns = 2:ncol(res))
  ) %>%
  cols_width(
    1 ~ px(165),
    2 ~ px(120),
    3 ~ px(50),
    4 ~ px(95),
    5 ~ px(120),
    6 ~ px(50),
    7 ~ px(95),
    8 ~ px(120),
    9 ~ px(50),
    10 ~ px(95),
    11 ~ px(120),
    12 ~ px(50),
    13 ~ px(95)
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
    source_notes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(
      rows = nrow(res)
    )
  ) %>%
  tab_source_note(
    source_note = "Note: [1] final model."
  ) %>%
  tab_source_note(
    source_note = "Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, SEM = Structural Equation Modeling, CI = Confidence Interval, ICU = Intensive Care Unit, HD = High Dependency."
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = which(rowSums(is.na(res[, 2:ncol(res)]) | res[, 2:ncol(res)] == "", na.rm = TRUE) == ncol(res)-1),
      columns = 1
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
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
      font = "Times New Roman",
      size = px(10)
    ),
    locations = cells_source_notes()
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",   
      size = px(10),             
      weight = "bold"           
    ),
    locations = cells_column_spanners(spanners = c("vap", "vap_sem_1", "vap_sem_2",
                                                   "bsi1", "bsi1_sem_1", "bsi1_sem_2"))
  ) %>%
  text_transform(
    locations = cells_body(
      columns = Characteristics,
      rows = Characteristics == "Acinetobacter spp."
    ),
    fn = function(x) html("<b><i>Acinetobacter</i> spp.</b>")
  )

print(table)

# Save
gtsave(data = table, filename = "output/table/table_analysis_sem.html")
