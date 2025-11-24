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
res3 <- readRDS("data/attfbis_table_multi_3_sub.RData")
#
res_vap <- cbind.data.frame(res1[[1]][,c(1,5,6,10)], 
                            res2[[1]][,c(5,6,10)],
                            res3[[1]][,c(5,6,10)])

res_bsi_1 <- cbind.data.frame(res1[[2]][,c(5,6,10)])

res_bsi_2 <- cbind.data.frame(res1[[3]][,c(5,6,10)], 
                              res2[[3]][,c(5,6,10)])

res <- cbind.data.frame(res_vap, res_bsi_1, res_bsi_2)

# Set col names
col_names <- as.character(unlist(res[1, ]))
res <- res[-1, ]
rownames(res) <- NULL

groups <- c("vap_1", "vap_2", "vap_3", 
            "bsi1_1", 
            "bsi2_1", "bsi2_2")

col_parts <- c("Unstandardized coefficient (95%CI)", "P value", "Standardized coefficient (95%CI)")

col_names <- c("Characteristics", unlist(lapply(groups, function(g) paste(col_parts, g, sep = "_"))))

colnames(res) <- col_names
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
    columns = 2:10,
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
    label = md("**<span style = 'color:black;'>SEM</span><sup> </sup>**"),
    columns = 5:7,
    level = 1,
    id = "vap_sem_2"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup>[1]</sup>**"),
    columns = 8:10,
    level = 1,
    id = "vap_sem_3"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Hospital-acquired BSI</span></sup>**"),
    columns = 11:13,
    level = 2,
    id = "bsi1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup>[1]</sup>**"),
    columns = 11:13,
    level = 1,
    id = "bsi1_sem_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Healthcare-associated BSI</span></sup>**"),
    columns = 14:19,
    level = 2,
    id = "bsi2"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup> </sup>**"),
    columns = 14:16,
    level = 1,
    id = "bsi2_sem_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>SEM</span><sup>[1]</sup>**"),
    columns = 17:19,
    level = 1,
    id = "bsi2_sem_2"
  ) %>%
  cols_label(
    `Unstandardized coefficient (95%CI)_vap_1` = md("Unstandardized coefficient (95%CI))"),
    `P value_vap_1` = md("*p*"),
    `Standardized coefficient (95%CI)_vap_1` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_vap_2` = md("Unstandardized coefficient (95%CI)"),
    `P value_vap_2` = md("*p*"),
    `Standardized coefficient (95%CI)_vap_2` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_vap_3` = md("Unstandardized coefficient (95%CI)"),
    `P value_vap_3` = md("*p*"),
    `Standardized coefficient (95%CI)_vap_3` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_bsi1_1` = md("Unstandardized coefficient (95%CI))"),
    `P value_bsi1_1` = md("*p*"),
    `Standardized coefficient (95%CI)_bsi1_1` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_bsi2_1` = md("Unstandardized coefficient (95%CI))"),
    `P value_bsi2_1` = md("*p*"),
    `Standardized coefficient (95%CI)_bsi2_1` = md("Standardized coefficient (95%CI)"),
    `Unstandardized coefficient (95%CI)_bsi2_2` = md("Unstandardized coefficient (95%CI)"),
    `P value_bsi2_2` = md("*p*"),
    `Standardized coefficient (95%CI)_bsi2_2` = md("Standardized coefficient (95%CI)")
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
    c(2, 5, 8, 11, 14, 17) ~ px(120),
    c(3, 6, 9, 12, 15, 18) ~ px(50),
    c(4, 7, 10, 13, 16, 19) ~ px(95)
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
    locations = cells_column_spanners(spanners = c("vap", "vap_sem_1", "vap_sem_2", "vap_sem_3",
                                                   "bsi1", "bsi1_sem_1",
                                                   "bsi2", "bsi2_sem_1", "bsi2_sem_2"))
  ) 

print(table)

#
gtsave(data = table, filename = "output/table/table_analysis_sem.html")
