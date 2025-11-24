# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 MASS)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_fbis.RData")

# Enterobacterales, carbapenem-resistant
df <- as.data.frame(df_all[[3]])

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

# Y as ordinal
df$fbis_score <- factor(df$fbis_score, ordered = TRUE)

# Model
formula <- fbis_score ~
  age_new + sex + country_region + country_income +
  hpd_admreason + severity_score_scale + icu_hd_ap + 
  pathogen_combined_types +
  ent_car  

# 
df_result <- bootstrap_prob_per_obs(
  df_input = df,
  formula  = formula,
  var_name = "ent_car",  
  n_boot   = 1000
)

# Mean by score_level
prob_summary <- function(df_subset, group_label) {
  
  df_subset %>%
    group_by(score_level) %>%
    summarise(
      mean_prob = mean(mean_prob, na.rm = TRUE),
      lower_ci  = mean(lower_ci, na.rm = TRUE),
      upper_ci  = mean(upper_ci, na.rm = TRUE)
    ) %>%
    mutate(
      fbis_score = as.numeric(score_level),
      group      = group_label
    ) %>%
    dplyr::select(fbis_score, mean_prob, lower_ci, upper_ci, group)
}


# 
df_plot <- rbind(
  prob_summary(
    df_result[grepl("susceptible$", df_result$ent_car), ], 
    "Carbapenem-susceptible"
  ),
  prob_summary(
    df_result[grepl("resistant$", df_result$ent_car), ], 
    "Carbapenem-resistant"
  )
)

# N total
N_group <- df_result %>%
  group_by(ent_car) %>%
  summarise(N = n_distinct(recordid)) %>%
  mutate(group = ifelse(
    grepl("susceptible$", ent_car),
    "Carbapenem-susceptible",
    "Carbapenem-resistant"
  )) %>%
  dplyr::select(group, N)

df_expected_count <- df_plot %>%
  left_join(N_group, by = "group") %>%
  mutate(
    expected_count = round(mean_prob * N, 0),
  ) %>%
  dplyr::select(fbis_score, expected_count, group)

print(df_expected_count)


# Save data
saveRDS(df_plot, "data/data_plot_3.RData")
saveRDS(df_expected_count, "data/df_count_3.RData")
###