# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 tidycmprsk, 
                 ggsurvfit,
                 ggplot2,
                 grid,
                 patchwork,
                 tidyverse,
                 survminer,
                 survival,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_data <- readRDS("data/rp_data.RData")

# 
df_data[[1]]$duration_vap_recur <- NULL
df_data[[2]]$duration_bsi_recur <- NULL

df <- rbind.data.frame(df_data[[1]], df_data[[2]]) %>%
  filter(time > 0)

# Load fonts
loadfonts()

# Convert outcome variable to numeric
df$event <- as.numeric(df$event)

#
df$time[df$time == 90 & df$event %in% c(1,2)] <- 89.999

# Set factor levels and labels for event
df$event <- factor(df$event,
                   levels = c(0, 2, 1),
                   labels = c("censor", "All-cause death", "Recurrence"))

# Calculate cumulative incidence using cuminc() which produces a tidycuminc object
cuminc <- cuminc(Surv(time, event) ~ infection_types, data = df)

# Convert the cuminc object into a tidy format using tidy_cuminc
cuminc_data <- tidy_cuminc(cuminc)

# Ensure that the data is sorted by time within each strata
cuminc_data <- cuminc_data %>%
  arrange(strata, time)

# Manually calculate n.risk (number at risk) by subtracting the cumulative number of events from the original number at risk
cuminc_data <- cuminc_data %>%
  group_by(strata) %>%
  mutate(n.risk = first(n.risk) - cumsum(n.event)) %>%  
  ungroup()

# Customize gg_cuminc plot
gg_cuminc <- cuminc |>
  ggcuminc(outcome = c("Recurrence", "All-cause death")) +
  scale_ggsurvfit() +
  scale_x_continuous(breaks = seq(0, 90, by = 15), limits = c(0, 90)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",  
    legend.box = "vertical",     
    legend.box.just = "center",
    legend.spacing.y = unit(0, "pt"),
    axis.text = element_text(size = 10, color = "black", family = "Times New Roman"),
    axis.title = element_text(size = 10, color = "black", family = "Times New Roman"),
    legend.text = element_text(size = 10, color = "black", family = "Times New Roman"),
    legend.title = element_text(size = 10, color = "black", face = "plain", family = "Times New Roman")
  ) +
  labs(x = "Follow-up duration since the first specimen date (days)",
       y = "Cumulative incidence") +
  add_confidence_interval() +
  scale_color_manual(
    name = "Infection syndromes",
    values = c("VAP" = "#1A75BB", 
               "Hospital-acquired BSI" = "#009344", 
               "Healthcare-associated BSI" = "#F05A28"),
    breaks = c("VAP", 
               "Hospital-acquired BSI", 
               "Healthcare-associated BSI")
  ) +
  scale_fill_manual(
    name = "Infection syndromes",
    values = c("VAP" = "#1A75BB", 
               "Hospital-acquired BSI" = "#009344", 
               "Healthcare-associated BSI" = "#F05A28"),
    breaks = c("VAP", 
               "Hospital-acquired BSI", 
               "Healthcare-associated BSI")
  ) +
  scale_linetype_manual(
    name = "Events",
    values = c("Recurrence" = "solid", "All-cause death" = "dashed"),
    breaks = c("Recurrence", "All-cause death"),
    labels = c("Recurrence", "All-cause death")
  ) +
  annotate("text", x = 8.2, y = 0.35, label = "p < 0.001 (all-cause death)",
           size = 3.5, color = "black", family = "Times New Roman") +
  annotate("text", x = 6.3, y = 0.40, label = "p < 0.001 (recurrence)", 
           size = 3.5, color = "black", family = "Times New Roman") +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    fill = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    linetype = guide_legend(nrow = 1, byrow = TRUE, order = 2)
  )

# Build the risk table plot with custom facet strip colors
gg_risktable <- 
  cuminc_data |>
  filter(time %in% seq(0, 90, by = 15)) |>
  select(outcome, strata, time, n.risk, cum.event) %>%
  {
    distinct(., strata, time, n.risk) |> 
      mutate(outcome = "At risk") |> 
      rename(stat = n.risk) |> 
      bind_rows(select(., outcome, strata, time, stat = cum.event))
  } |>
  mutate(outcome = recode(outcome, 
                          "All-cause death" = "Adeath", 
                          "Recurrence" = "Recurr", 
                          "At risk" = "At risk")) |> 
  mutate(outcome = factor(outcome, levels = c("Adeath", "Recurr", "At risk"))) |>
  # Set the strata factor levels in the correct order
  mutate(strata = factor(strata, levels = c("VAP", 
                                            "Hospital-acquired BSI", 
                                            "Healthcare-associated BSI"))) |>
  ggplot(aes(x = time, y = factor(outcome), label = stat, linetype = outcome)) +  
  geom_text(size = 3.5, color = "black", family = "Times New Roman",
            position = position_dodge(width = 0.9), 
            check_overlap = TRUE) +
  labs(y = NULL, x = NULL) +
  facet_grid(strata ~ .) +
  coord_cartesian(xlim = c(0, 90)) +
  scale_x_continuous(breaks = seq(0, 90, by = 15)) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 0),
    strip.placement = "outside",
    axis.text.y = element_text(size = 10, color = "black", family = "Times New Roman"),
    axis.text.x = element_blank()
  ) +
  scale_color_manual(
    name = "Infection syndromes",
    values = c("VAP" = "#1A75BB", 
               "Hospital-acquired BSI" = "#009344", 
               "Healthcare-associated BSI" = "#F05A28"),
    breaks = c("VAP", 
               "Hospital-acquired BSI", 
               "Healthcare-associated BSI")
  ) +
  scale_linetype_manual(
    name = "Events", 
    values = c("Adeath" = "dashed", "Recur" = "solid", "At risk" = "solid"),  
    breaks = c("Adeath", "Recurr", "At risk"),  
    labels = c("All-cause death", "Recurrence", "At risk")
  )

# Convert the ggplot object to a gtable
g <- ggplot_gtable(ggplot_build(gg_risktable))

# Identify the facet strips
stripr <- which(grepl('strip-r', g$layout$name))

# Define the colors for the facet strips
fills <- c(adjustcolor("#1A75BB", alpha.f = 0.7), 
           adjustcolor("#009344", alpha.f = 0.7),
           adjustcolor("#F05A28", alpha.f = 0.7))

# Modify the background color of the facet strips
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k + 1
}

# Align and combine the plots
gg_combined <- list(gg_cuminc, g)

# Visualize the combined plots
combined_plot <- wrap_plots(
  gg_combined[[1]], 
  gg_combined[[2]], 
  ncol = 1,
  heights = c(1, 0.8)
) + 
  plot_layout(widths = c(1.5)) + 
  plot_annotation(
    theme = theme(plot.margin = margin(10, 5, 5, 50))
  )

# Print the combined plot
print(combined_plot)

# Save the figure to a PDF file
cairo_pdf(file = "output/Recurrence_CIF.pdf", 
          width = 8, height = 5.5)
print(combined_plot)
dev.off()
###