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
res1 <- readRDS("data/attfbis_value_uni.RData")
res2 <- readRDS("data/attfbis_value_multi_1.RData")
res3 <- readRDS("data/attfbis_value_multi_2.RData")
#
res_vap <- cbind.data.frame(res1[[1]], 
                            res2[[1]][,c(2,3)],
                            res3[[1]][,c(2,3)])

res_bsi_1 <- cbind.data.frame(res1[[2]][,c(2,3)], 
                              res2[[2]][,c(2,3)])

res_bsi_2 <- cbind.data.frame(res1[[3]][,c(2,3)], 
                              res2[[3]][,c(2,3)],
                              res3[[3]][,c(2,3)])

res <- cbind.data.frame(res_vap, res_bsi_1, res_bsi_2)

# Set col names
col_names <- as.character(unlist(res[1, ]))
res <- res[-1, ]
rownames(res) <- NULL

colnames(res) <- c("Characteristics", 
                   "Crude coefficient (95%CI)_uni_vap", "P value_uni_vap",
                   "Adjusted coefficient (95%CI)_multi_vap_1", "P value_multi_vap_1",
                   "Adjusted coefficient (95%CI)_multi_vap_2", "P value_multi_vap_2",
                   "Crude coefficient (95%CI)_uni_bsi1", "P value_uni_bsi1",
                   "Adjusted coefficient (95%CI)_multi_bsi1_1", "P value_multi_bsi1_1",
                   "Crude coefficient (95%CI)_uni_bsi2", "P value_uni_bsi2",
                   "Adjusted coefficient (95%CI)_multi_bsi2_1", "P value_multi_bsi2_1",
                   "Adjusted coefficient (95%CI)_multi_bsi2_2", "P value_multi_bsi2_2")
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
  mutate(across(everything(), ~ ifelse(.x == "NA", "", .x))) %>%
  mutate(across(everything(), as.character)) %>%
  pmap_dfr(process_missing_values)

# gt table
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
    label = md("**<span style = 'color:black;'>Univariable analysis</span><sup> </sup>**"),
    columns = 2:3,
    level = 1,
    id = "uni_vap"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Multivariable analysis</span><sup> </sup>**"),
    columns = 4:5,
    level = 1,
    id = "multi_vap_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Multivariable analysis</span><sup>[1]</sup>**"),
    columns = 6:7,
    level = 1,
    id = "multi_vap_2"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Hospital-acquired BSI</span></sup>**"),
    columns = 8:11,
    level = 2,
    id = "bsi1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Univariable analysis</span><sup> </sup>**"),
    columns = 8:9,
    level = 1,
    id = "uni_bsi1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Multivariable analysis</span><sup>[1]</sup>**"),
    columns = 10:11,
    level = 1,
    id = "multi_bsi1_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Healthcare-associated BSI</span></sup>**"),
    columns = 12:ncol(res),
    level = 2,
    id = "bsi2"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Univariable analysis</span><sup> </sup>**"),
    columns = 12:13,
    level = 1,
    id = "uni_bsi2"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Multivariable analysis</span><sup> </sup>**"),
    columns = 14:15,
    level = 1,
    id = "multi_bsi2_1"
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Multivariable analysis</span><sup>[1]</sup>**"),
    columns = 16:17,
    level = 1,
    id = "multi_bsi2_2"
  ) %>%
  tab_style(
    style = cell_text(align = "center", v_align = "middle"),
    locations = cells_column_labels(
      columns = c(
        `P value_uni_vap`,
        `P value_multi_vap_1`,
        `P value_multi_vap_2`,
        `P value_uni_bsi1`,
        `P value_multi_bsi1_1`,
        `P value_uni_bsi2`,
        `P value_multi_bsi2_1`,
        `P value_multi_bsi2_2`
      )
    )
  ) %>%
  cols_label(
    `Crude coefficient (95%CI)_uni_vap` = md("Crude coefficient (95%CI)"),
    `P value_uni_vap` = md("*p*"),
    `Adjusted coefficient (95%CI)_multi_vap_1` = md("Adjusted coefficient (95%CI)"),
    `P value_multi_vap_1` = md("*p*"),
    `Adjusted coefficient (95%CI)_multi_vap_2` = md("Adjusted coefficient (95%CI)"),
    `P value_multi_vap_2` = md("*p*"),
    `Crude coefficient (95%CI)_uni_bsi1` = md("Crude coefficient (95%CI)"),
    `P value_uni_bsi1` = md("*p*"),
    `Adjusted coefficient (95%CI)_multi_bsi1_1` = md("Adjusted coefficient (95%CI)"),
    `P value_multi_bsi1_1` = md("*p*"),
    `Crude coefficient (95%CI)_uni_bsi2` = md("Crude coefficient (95%CI)"),
    `P value_uni_bsi2` = md("*p*"),
    `Adjusted coefficient (95%CI)_multi_bsi2_1` = md("Adjusted coefficient (95%CI)"),
    `P value_multi_bsi2_1` = md("*p*"),
    `Adjusted coefficient (95%CI)_multi_bsi2_2` = md("Adjusted coefficient (95%CI)"),
    `P value_multi_bsi2_2` = md("*p*")
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
  cols_width(
    1 ~ px(200),
    c(seq(2, 16, 2)) ~ px(120),
    c(seq(3, 17, 2)) ~ px(60)
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
    source_note = "Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, CI = Confidence Interval, ICU = Intensive Care Unit, HD = High Dependency."
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
    locations = cells_column_spanners(spanners = c("vap",
                                                   "uni_vap", 
                                                   "multi_vap_1", "multi_vap_2",
                                                   "bsi1",
                                                   "uni_bsi1",
                                                   "multi_bsi1_1",
                                                   "bsi2",
                                                   "uni_bsi2",
                                                   "multi_bsi2_1",
                                                   "multi_bsi2_2"))
  ) 

print(table)

# Save
gtsave(data = table, filename = "output/table/table_analysis_value.html")
