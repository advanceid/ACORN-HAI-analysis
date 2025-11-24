# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 survival,
                 survminer,
                 ggplot2,
                 gridExtra,
                 ggpubr,
                 extrafont,
                 patchwork,
                 Cairo,
                 stringr,
                 grid)
})

# Set working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_first28.RData")

# Load fonts
loadfonts()

# Variables
var <- c("aci_car", "ent_thir", "ent_car", 
         "pse_car", "entc_van", "sa_meth", 
         "mdr_gnb", "mdr", "amr")

# Raw labels for group assignment
raw_labels <- list(
  c("CSA", "CRA"),
  c("3GCSE", "3GCRE"),
  c("CSE", "CRE"),
  c("CSP", "CRP"),
  c("VSE", "VRE"),
  c("MSSA", "MRSA"),
  c("Non-MDR-GNB", "MDR-GNB"),
  c("Non-MDR", "MDR"),
  c("Non-AMR", "AMR")
)

# Create padded labels for visual alignment
max_len <- max(nchar(unlist(raw_labels)))
labels <- lapply(raw_labels, function(x) str_pad(x, max_len, side = "left"))

# Initialize list to store plots
combined_plot <- list()

# Loop over each dataset
for (i in 1:6) {
  df <- df_all[[i]]
  ris_name <- var[i]
  
  # Define groupings
  df$ris_group <- factor(
    ifelse(grepl("resistant$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][2],
           ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][1], NA)),
    levels = raw_labels[[i]]
  )
  
  df$ris_label <- factor(
    ifelse(grepl("resistant$", df[[ris_name]], ignore.case = TRUE), labels[[i]][2],
           ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), labels[[i]][1], NA)),
    levels = labels[[i]]
  )
  
  df$time[df$time == 28 & df$event == 1] <- 27.999
  
  # Fit KM model
  fit <- survfit(Surv(time, event) ~ ris_group, data = df)
  
  # p-value
  p_model <- survdiff(Surv(time, event) ~ ris_group, data = df)
  p_value <- 1 - pchisq(p_model$chisq, df = length(p_model$n) - 1)
  p_label <- ifelse(p_value < 0.001, "< 0.001",
                    ifelse(p_value < 0.01, paste0("= ", sprintf("%.3f", p_value)),
                           paste0("= ", sprintf("%.2f", p_value))))
  
  # per-plot margins for y text/title 
  text_r  <- c(13, 15, 10, 10, 10, 10) 
  title_r <- c(15, 20, 15, 15, 15, 15) 
  text_r_table <- c(14, 12, 12, 12, 12, 12) 
  
  # Plot survival
  plt_fct <- ggsurvplot(
    fit,
    data = df,
    pval = FALSE,
    palette = c("#EFC000FF", "#0073c2FF"),
    conf.int = TRUE,
    linetype = 1,
    ggtheme = theme_minimal(base_family = "Times New Roman") +
      theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 16, family = "Times New Roman", hjust = 1.2),
        axis.text.y = element_text(color = "black", size = 16, family = "Times New Roman", 
                                   margin = margin(r = text_r[i])),
        axis.title.y = element_text(color = "black", size = 16, family = "Times New Roman", 
                                    margin = margin(r = title_r[i])),
        axis.title.x = element_text(color = "black", size = 16, family = "Times New Roman", margin = margin(t = 10)),
        panel.spacing = unit(0.6, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        strip.text.x = element_text(size = 16, margin = margin(t = 6, b = 6)),
        strip.background.x = element_rect(color = "black", fill = "gray90", linewidth = NA)
      ),
    title = paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1]),
    xlab = "Follow-up time since infection onset (days)",
    ylab = "Survival probability",
    font.tickslab = c(16, "plain"),
    font.x = c(16, "plain"),
    font.y = c(16, "plain"),
    break.time.by = 7,
    xlim = c(0, 28),
    legend = "none"
  )
  
  # Risk table
  ggsurv <- ggsurvplot(fit, df, risk.table = TRUE, break.time.by = 7, xlim = c(0, 28), tables.theme = theme_minimal())
  
  tbl_data <- ggsurv$table$data %>%
    mutate(
      ris_group = gsub("ris_group=", "", strata),
      ris_group = factor(ris_group, levels = raw_labels[[i]])
    ) %>%
    filter(!is.na(n.risk), time <= 28)
  
  tbl_fct <- ggplot(tbl_data, aes(time, ris_group, color = ris_group)) +
    geom_text(aes(label = n.risk), size = 6, family = "Times New Roman") +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
      strip.background = element_rect(fill = "white", linewidth = 1),
      axis.text.y = element_text(color = "black", size = 16, family = "Times New Roman", 
                                 margin = margin(r = text_r_table[i])),
      axis.text.x = element_blank(),
      strip.text = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.6, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")
    ) +
    scale_y_discrete(labels = labels[[i]]) +
    scale_x_continuous(
      limits = c(0, 28),
      breaks = c(0, 7, 14, 21, 28)
    ) +
    scale_color_manual(
      values = setNames(c("#EFC000FF", "#0073c2FF"), raw_labels[[i]]),
      labels = labels[[i]]
    ) +
    ggtitle(" ") + xlab("") + ylab("") + guides(color = "none")
  
  # Add title and p-value
  plt_fct$plot <- plt_fct$plot + ggtitle(paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1])) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, family = "Times New Roman")) +
    annotate("text", x = 5, y = 0.15,
             label = paste0("p ", p_label),
             size = 5.2, hjust = 1,
             family = "Times New Roman", fontface = "plain")
  
  # Combine plot + risk table
  grob_combined <- arrangeGrob(plt_fct$plot, tbl_fct, heights = c(7.5, 2.5))
  
  grob_padded <- gTree(
    children = gList(grob_combined),
    vp = viewport(layout = grid.layout(1, 1),
                  width = unit(1, "npc") - unit(1, "cm"),
                  height = unit(1, "npc") - unit(0.5, "cm"))
  )
  
  # Store the final plot in the list
  combined_plot[[i]] <- grob_padded

}

# Save
for (i in 1:6) {
  CairoPDF(file = paste0("output/pdf/att28_death_KM_0", i, "_overall.pdf"), width = 7, height = 5)
  grid.draw(combined_plot[[i]])
  dev.off()
}

# =========
# Helper: compute crude mortality by group
# =========
compute_crude <- function(df_all, var, raw_labels) {
  stopifnot(length(df_all) == length(var), length(var) == length(raw_labels))
  
  out_list <- vector("list", length(var))
  
  for (i in seq_along(var)) {
    df <- df_all[[i]]
    ris_name <- var[i]
    
    # Grouping consistent with your KM plots
    df$ris_group <- ifelse(
      grepl("resistant$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][2],
      ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][1], NA)
    )
    df$ris_group <- factor(df$ris_group, levels = raw_labels[[i]])
    
    # Summarize by ris_group
    tmp <- df %>%
      dplyr::filter(!is.na(ris_group)) %>%
      dplyr::group_by(ris_group) %>%
      dplyr::summarise(
        total   = dplyr::n(),
        deaths  = sum(event == 1, na.rm = TRUE),   # 28-day deaths
        crude   = deaths / total,
        crude_pct = scales::percent(crude, accuracy = 0.1),
        .groups = "drop"
      ) %>%
      dplyr::mutate(pattern = paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1]))
    
    out_list[[i]] <- tmp
  }
  
  dplyr::bind_rows(out_list) %>%
    dplyr::select(pattern, ris_group, total, deaths, crude, crude_pct)
}

crude_mortality <- compute_crude(df_all, var, raw_labels)
crude_mortality


