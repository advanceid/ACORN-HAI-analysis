# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 MASS,
                 ggplot2,
                 scales,
                 Cairo,
                 extrafont,
                 patchwork)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_fbis.RData")

# Load fonts
loadfonts()

# Enterobacterales, third-generation cephalosporin-resistant
df_all <- as.data.frame(df_all[[2]])

# Split by infection types
df_used <- split(df_all, df_all$infection_types)

# CI
bootstrap_prob_per_obs <- function(df_input, formula, var_name, n_boot = 500) {

  if(!is.ordered(df_input$fbis_score)){
    df_input$fbis_score <- factor(df_input$fbis_score, ordered = TRUE)
  }
  
  n_obs <- nrow(df_input)
  score_levels <- levels(df_input$fbis_score)
  K <- length(score_levels)
  
  prob_array <- array(
    NA,
    dim = c(n_obs, K, n_boot),
    dimnames = list(NULL, score_levels, NULL)
  )
  
  # Model
  for (b in seq_len(n_boot)) {
    
    boot_sample <- df_input[sample(n_obs, replace = TRUE), ]
    
    if (length(unique(boot_sample$fbis_score)) < K) next
    
    boot_model <- try(
      MASS::polr(formula, data = boot_sample, method = "logistic", Hess = TRUE),
      silent = TRUE
    )
    if (inherits(boot_model, "try-error")) next
    
    # 
    pred_probs <- try(
      predict(boot_model, newdata = df_input, type = "probs"),
      silent = TRUE
    )
    if (inherits(pred_probs, "try-error")) next
    
    # 
    prob_array[, , b] <- pred_probs
  }
  
  # 
  valid_boots <- which(!is.na(prob_array[1, 1, ]))
  cat("Successful bootstraps:", length(valid_boots), "of", n_boot, "\n")
  
  # Save
  results_list <- vector("list", K)
  for (k in seq_len(K)) {
    probs_k <- prob_array[, k, valid_boots, drop = FALSE]
    
    df_k <- data.frame(
      recordid    = df_input$recordid,
      obs_id      = seq_len(n_obs),
      score_level = score_levels[k],
      mean_prob   = rowMeans(probs_k, na.rm = TRUE),
      lower_ci    = apply(probs_k, 1, quantile, probs = 0.025, na.rm = TRUE),
      upper_ci    = apply(probs_k, 1, quantile, probs = 0.975, na.rm = TRUE)
    )
    
    # 
    df_k[[var_name]] <- df_input[[var_name]]
    
    results_list[[k]] <- df_k
  }
  
  final_result <- do.call(rbind, results_list)
  rownames(final_result) <- NULL
  return(final_result)
}

# Result
set.seed(123)
df_plot_list <- list()

for (i in 1:3) {
  df <- df_used[[i]]
  
  # Save infection type information
  infection_type <- unique(df$infection_types)
  
  # Ensure fbis_score is an ordered factor
  df$fbis_score <- factor(df$fbis_score, ordered = TRUE)
  
  # Define model formula based on infection type
  if (i == 1) {
    formula <- fbis_score ~ 
      age_new + sex + 
      country_income + hpd_admreason + 
      severity_score_scale + icu_hd_ap +
      ent_thir  
  } else if (i == 2) {
    formula <- fbis_score ~ 
      age_new + sex + 
      country_income +
      comorbidities_CCI +
      severity_score_scale + icu_hd_ap + 
      ent_thir  
  } else if (i == 3) {
    formula <- fbis_score ~ 
      age_new + sex + comorbidities_CCI +
      icu_hd_ap + 
      pathogen_combined_types +
      ent_thir 
  }
  
  # Define a summary function that now also adds the infection_type column
  prob_summary <- function(df_subset, group_label, infection_type) {
    df_subset %>%
      group_by(score_level) %>%
      summarise(
        mean_prob = mean(mean_prob, na.rm = TRUE),
        lower_ci = mean(lower_ci, na.rm = TRUE),
        upper_ci = mean(upper_ci, na.rm = TRUE)
      ) %>%
      mutate(
        fbis_score = as.numeric(score_level),
        group = group_label,
        infection_type = infection_type
      ) %>%
      dplyr::select(infection_type, fbis_score, mean_prob, lower_ci, upper_ci, group)
  }
  
  # Perform bootstrap
  df_result <- bootstrap_prob_per_obs(
    df_input = df,
    formula  = formula,
    var_name = "ent_thir",  
    n_boot   = 1000
  )
  
  # For each group, filter based on ent_thir values and summarize
  summary_S <- prob_summary(
    df_result[grepl("susceptible$", df_result$ent_thir), ], 
    "Third-generation cephalosporin-susceptible", infection_type
  )
  summary_R <- prob_summary(
    df_result[grepl("resistant$", df_result$ent_thir), ], 
    "Third-generation cephalosporin-resistant", infection_type
  )
  
  df_plot_list[[i]] <- rbind(summary_S, summary_R)
}

# Combine for plot
df_plot <- do.call(rbind, df_plot_list)

# Plot with facets for each infection type
p <- ggplot(df_plot, aes(x = factor(fbis_score), 
                         y = mean_prob * 100, 
                         color = group, group = group)) +
  geom_point(size = 1.5) +  
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower_ci * 100, ymax = upper_ci * 100),  
                width = 0.4, linewidth = 0.5) + 
  scale_y_continuous(
    labels = label_percent(scale = 1),  
    breaks = seq(0, 100, by = 10), 
  ) +
  labs(
    title = "",  
    x = "FBIS score",
    y = "Predicted probability"
  ) +
  facet_wrap(~ infection_type) +  
  scale_color_manual(
    name = "Enterobacterales",  
    values = c("#0073c2FF", "#EFC000FF"),  
    labels = c("Third-generation cephalosporin-resistant", 
               "Third-generation cephalosporin-susceptible")  
  )  +  
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 10), 
    strip.text = element_text(family = "Times New Roman", size = 10),
    axis.title.x = element_text(family = "Times New Roman", size = 10, 
                                margin = margin(t = 5)),
    axis.title.y = element_text(family = "Times New Roman", size = 10, 
                                margin = margin(r = 5)),
    axis.text.y = element_text(size = 10, family = "Times New Roman", 
                               color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(size = 10, family = "Times New Roman", 
                               color = "black", margin = margin(t = 3)),
    legend.title = element_text(size = 10, 
                                family = "Times New Roman"), 
    legend.text = element_text(size = 10, family = "Times New Roman"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.3, color = "black", fill = NA),
    strip.background = element_rect(linewidth = 0.3, color = "black", fill = "gray80"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(t = -7, r = 0, b = 0, l = 0),
    plot.margin = margin(t = -10, r = 0, b = 0, l = 0)
  ) 

p

# Save
CairoPDF("output/pdf/fbis_2_sub.pdf", width = 9, height = 4)
print(p)
dev.off()
