# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 ggplot2,
                 scales,
                 patchwork,
                 Cairo,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df1 <- readRDS("1_car_aci/data/data_plot_1.RData")
df2 <- readRDS("2_thir_ent/data/data_plot_2.RData")
df3 <- readRDS("3_car_ent/data/data_plot_3.RData")

df_count_1 <- readRDS("1_car_aci/data/df_count_1.RData")
df_count_2 <- readRDS("2_thir_ent/data/df_count_2.RData")
df_count_3 <- readRDS("3_car_ent/data/df_count_3.RData")
# Load fonts
loadfonts()

# Counts
df_count_1$group <- ifelse(grepl("susceptible$", df_count_1$group),
                           "CSA",
                           ifelse(grepl("resistant$", df_count_1$group),
                                  "CRA",
                                  df_count_1$group))


df_count_2$group <- ifelse(grepl("susceptible$", df_count_2$group),
                           "3GCSE",
                           ifelse(grepl("resistant$", df_count_2$group),
                                  "3GCRE",
                                  df_count_2$group))


df_count_3$group <- ifelse(grepl("susceptible$", df_count_3$group),
                           "CSE",
                           ifelse(grepl("resistant$", df_count_3$group),
                                  "CRE",
                                  df_count_3$group))

df_counts <- rbind(df_count_1, df_count_2, df_count_3)


# Merge
merge_counts <- function(df_pred, df_counts, groups_raw, groups_named) {
  
  group_map <- setNames(groups_named, groups_raw)
  
  df_counts_sub <- df_counts %>%
    filter(group %in% groups_raw) %>%
    mutate(
      group = recode(group, !!!group_map),
      fbis_score = factor(fbis_score, levels = 1:7, ordered = TRUE)
    )
  
  df_pred %>%
    mutate(
      fbis_score = factor(fbis_score, levels = 1:7, ordered = TRUE)
    ) %>%
    left_join(df_counts_sub, by = c("fbis_score", "group"))
}
# df1
df1 <- merge_counts(
  df1, df_counts,
  groups_raw = c("CRA", "CSA"),
  groups_named = c("Carbapenem-resistant", "Carbapenem-susceptible")
)

# df2
df2 <- merge_counts(
  df2, df_counts,
  groups_raw = c("3GCRE", "3GCSE"),
  groups_named = c("Third-generation cephalosporin-resistant", "Third-generation cephalosporin-susceptible")
)

# df3 
df3 <- merge_counts(
  df3, df_counts,
  groups_raw = c("CRE", "CSE"),
  groups_named = c("Carbapenem-resistant", "Carbapenem-susceptible")
)


# For plot
p_list <- list(df1, df2, df3)

p <- list()

system.time({
  for (i in 1:3) {
    df <- p_list[[i]]
    
    df$fbis_score <- as.factor(df$fbis_score)
    
    if (i == 1) {
      p_title <- "CRA versus CSA"
    } else if (i == 2) {
      p_title <- "3GCRE versus 3GCSE"
    } else if (i == 3) {
      p_title <- "CRE versus CSE"
    }
    
    df_p <- ggplot(df, aes(x = factor(fbis_score), 
                           y = mean_prob * 100, 
                           color = group, group = group)) +
      geom_point(size = 1.5) +  
      geom_line(linewidth = 0.5) +
      geom_errorbar(aes(ymin = lower_ci * 100, 
                        ymax = upper_ci * 100),  
                    width = 0.3, linewidth = 0.5) + 
      scale_y_continuous(
        breaks = seq(0, 100, by = 10),
        labels = function(x) ifelse(x %% 20 == 0, scales::label_percent(scale = 1)(x), ""),
        limits = c(0, 100)
      ) +
      labs(
        title = p_title,  
        x = " ",
        y = if (i == 2) "Predicted probability" else "" 
      ) +  
      scale_color_manual(
        name = " ",  
        values = c("#0073c2FF", "#EFC000FF"),  
        labels = c("Resistant (R)", "Susceptible (S)")  
      ) +  
      geom_text(
        data = df %>% filter(grepl("susceptible$", group)),
        aes(x = fbis_score, y = 85, label = paste0("S = ", expected_count)),
        color = "#EFC000FF",
        size = 3.5,
        family = "Times New Roman",
        hjust = 0,
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df %>% filter(grepl("resistant$", group)),
        aes(x = fbis_score, y = 68, label = paste0("R = ", expected_count)),
        color = "#0073c2FF",
        size = 3.5,
        family = "Times New Roman",
        hjust = 0,
        inherit.aes = FALSE
      ) +
      theme_bw() +
      theme(
        text = element_text(family = "Times New Roman", size = 12), 
        strip.text = element_text(family = "Times New Roman", size = 12),
        axis.text.y = if (i == 1) element_text(size = 12, family = "Times New Roman") else element_blank(),
        axis.text.x = element_text(size = 12, family = "Times New Roman", 
                                   color = "black", margin = margin(t = 3)),
        axis.title.x = element_text(family = "Times New Roman", size = 12, 
                                    margin = margin(t = 5)),
        plot.title = element_text(family = "Times New Roman",  
                                  size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, family = "Times New Roman"), 
        legend.text = element_text(size = 12, family = "Times New Roman"),
        legend.position = if (i == 2) "bottom" else "none",
        legend.direction = "horizontal",
        legend.margin = margin(t = -7, r = 0, b = 0, l = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
      ) +
      coord_flip() 
    
    p[[i]] <- df_p
    
  }
})


# FBIS labels
fbis_labels <- data.frame(
  fbis_score = factor(1:7, ordered = TRUE),
  severity = 1:7,  
  label = c(
    "On palliative care in terminal phases.",
    "Accommodated in a long-term ventilator unit.",
    "Hospitalized in an intensive care unit.",
    "Hospitalized but not requiring\nan intensive care unit.",
    "Out of hospital; significant disability;\nrequires assistance.",
    "Out of hospital; moderate signs or symptoms;\nunable to complete daily activities.",
    "Out of hospital; basically healthy;\nable to complete activities."
  )
)

# Same factor
df1$fbis_score <- as.factor(df1$fbis_score)
fbis_labels$fbis_score <- factor(fbis_labels$fbis_score, levels = levels(df1$fbis_score))

#
p_heatmap <- ggplot(fbis_labels, aes(x = 1, y = fbis_score, fill = severity)) +
  geom_tile(color = "white", linewidth = 0.5) + 
  geom_text(aes(label = label), color = "white", 
            fontface = "bold", size = 3.5, hjust = 0.5) + 
  scale_fill_gradientn(colors = viridis_pal(option = "plasma", end = 0.9, direction = 1)(7)) +  
  labs(
    title = "FBIS score",
    x = NULL, 
    y = NULL
  ) +
  scale_y_discrete(position = "right") + 
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(6, "mm"),
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", 
                              size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_heatmap


#
final_plot <- p_heatmap + p[[1]] + p[[2]] + p [[3]] + plot_layout(widths = c(2,2,2,2))
final_plot

# Save
CairoPDF("pdf/fbis_overall.pdf", width = 16, height = 7)
print(final_plot)
dev.off()
###