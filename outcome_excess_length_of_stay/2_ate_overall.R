# Clear environment
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 cobalt,
                 mets,
                 extrafont,
                 scales)
})

# Load fonts
loadfonts()

# Read in the list of data frames for different time cut-offs
all_sub <- readRDS("data/los_lists.RData")

# Resistance pattern names
pattern_names <- c(
  "CRA vs CSA", "3GCRE vs 3GCSE", "CRE vs CSE",
  "CRP vs CSP", "VRE vs VSE", "MRSA vs MSSA"
)

# Time cut-off identifiers
cutoffs <- c("t14", "t21", "t28", "t60", "t90")

# Prepare list to hold all results
all_results <- list()

# Candidate covariates for Propensity Score Model
candidate_vars_full <- c(
  "age_new",
  "sex",
  "comorbidities_CCI",
  "severity_score_scale",
  "icu_hd_ap",
  "pathogen_combined_types"
)

# Stepwise AIC selection function (robust version)
select_vars_logit <- function(data, outcome, candidates) {
  data <- data %>% filter(!is.na(.data[[outcome]]))
  if (nrow(data) < 50) return(character(0))
  
  # Filter candidates to only include those with variance
  valid_candidates <- candidates[sapply(candidates, function(var) {
    if (!(var %in% colnames(data))) return(FALSE)
    n_levels <- length(unique(data[[var]][!is.na(data[[var]])]))
    n_levels >= 2  
  })]
  
  if (length(valid_candidates) == 0) return(character(0))
  
  full_formula <- as.formula(paste(outcome, "~", paste(valid_candidates, collapse = " + ")))
  null_formula <- as.formula(paste(outcome, "~ 1"))
  
  model_full <- tryCatch({
    glm(full_formula, data = data, family = binomial)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(model_full)) return(character(0))
  
  model_step <- step(model_full,
                     scope = list(lower = null_formula, upper = full_formula),
                     direction = "both", trace = 0)
  
  vars <- attr(terms(model_step), "term.labels")
  return(vars)
}


# Main analysis loop for each cut-off
for (tc in cutoffs) {
  
  df_list <- all_sub[[tc]]
  t_int <- as.numeric(sub("t", "", tc))
  
  # List to hold results for the current cutoff
  tc_results_list <- list() 
  
  for (i in seq_along(df_list)) {
    df <- as.data.frame(df_list[[i]])
    pat_name <- pattern_names[i]
    
    n_total <- nrow(df)
    n_levels <- length(unique(df$ast_ris))
    tab_event <- table(df$event_new)
    n_evt1 <- ifelse("1" %in% names(tab_event), tab_event["1"], 0)
    
    n_analyzed <- n_total 
    
    if (n_total >= 100 && n_levels == 2 && n_evt1 >= 30) {
      
      set.seed(42)
      selected_vars <- tryCatch({
        select_vars_logit(df, "ast_ris", candidate_vars_full)
      }, error = function(e) character(0))
      
      required_vars <- intersect(c("age_new", "sex"), colnames(df))
      selected_vars <- union(required_vars, selected_vars)
      selected_vars <- selected_vars[selected_vars != ""]
      
      # 1. Define Propensity Score Model Formula
      if (length(selected_vars) == 0) {
        treat_formula <- as.formula("ast_ris ~ 1")
      } else {
        treat_formula <- as.formula(paste("ast_ris ~", paste(selected_vars, collapse = " + ")))
      }
      
      # PS Calculation
      df_to_use <- df 
      
      ps_model <- tryCatch({
        glm(treat_formula, data = df, family = binomial)
      }, error = function(e) {
        message(paste("PS Model Error:", e$message, "Skipping Capping for", pat_name, tc))
        return(NULL)
      })
      
      if (!is.null(ps_model)) {
        # Calculate PS
        ps_scores <- predict(ps_model, type = "response", na.action = na.pass)
        df$ps <- ps_scores
        
        # Define Capping Thresholds
        trim_lower <- quantile(df$ps, 0.025, na.rm = TRUE)
        trim_upper <- quantile(df$ps, 0.975, na.rm = TRUE)
        
        # Apply Capping/Winsorizing to the PS score: replace extreme values
        df$ps_capped <- pmin(pmax(df$ps, trim_lower), trim_upper)
        
        # Count how many patients had their PS value capped (for messaging)
        n_capped_count <- sum(df$ps != df$ps_capped, na.rm = TRUE)
        
        df_to_use <- df
        n_analyzed <- n_total 
        
        message(paste0("Capping applied (", tc, ", ", pat_name, 
                       "). PS values capped for: ", n_capped_count, 
                       " patients. Total N: ", n_analyzed, "."))
      }
      
      # Define Censoring model (based on selected variables)
      if (all(df_to_use$event_new != 0)) {
        cens_formula <- ~1
      } else {
        # Censoring model uses ast_ris and dynamically selected variables
        cens_vars <- paste(c("ast_ris", selected_vars), collapse = " + ")
        cens_formula <- as.formula(paste("~", cens_vars))
      }
      
      # Run RMST model on the full data (df_to_use)
      res1 <- tryCatch({
        fit <- resmeanATE(
          formula = Event(time_new, event_new) ~ ast_ris,
          data = df_to_use,
          cause = 1,
          time = t_int,
          model = "exp",
          treat.model = treat_formula,  
          augmentation = TRUE,
          cens.model = cens_formula,
          type = "II",
          prop.haz = FALSE, # Use non-PH to increase stability
          kaplan.meier = TRUE,
          se = TRUE
        )
        
        s <- summary(fit)$ateDR["treat:1-0", ]
        
        if (is.null(s) || all(is.na(s))) {
          message(paste("ATE estimate failed post-capping for:", pat_name, tc))
          return(NULL)
        }
        
        data.frame(
          Cutoff = tc,
          Pattern = pat_name,
          Total_N = n_total,
          N_analyzed = n_analyzed,
          Estimate = -round(as.numeric(s["Estimate"]), 0),
          CI_lower = -round(as.numeric(s["97.5%"]), 0),
          CI_upper = -round(as.numeric(s["2.5%"]), 0),
          P_value = as.numeric(s["P-value"]),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        message(paste("Fatal Error in resmeanATE (Post-Capping) for:", pat_name, tc, "-", e$message))
        return(NULL)
      })
      
      if (!is.null(res1)) {
        tc_results_list[[length(tc_results_list) + 1]] <- res1
      }
    }
  }
  
  # Combine and format results for the current cutoff
  tc_results <- bind_rows(tc_results_list)
  
  if (!is.null(tc_results) && nrow(tc_results) > 0) {
    tc_results <- tc_results %>%
      mutate(
        Group = paste0("Total = ", Total_N), 
        Estimate_CI = sprintf("%.0f [%.0f, %.0f]", Estimate, CI_lower, CI_upper),
        P_value = case_when(
          P_value < 0.001 ~ "<0.001",
          P_value < 0.01  ~ sprintf("%.3f", P_value),
          TRUE ~ sprintf("%.2f", P_value)
        )
      ) %>%
      select(Cutoff, Pattern, Group, Estimate_CI, P_value)
    
    all_results[[tc]] <- tc_results
  }
}

# Combine all results and export
excess_los <- dplyr::bind_rows(all_results)
excess_los$group_by <- "Overall"

# Save
write.csv(excess_los, "output/excess_los_overall.csv", row.names = FALSE)
saveRDS(excess_los, "data/excess_los_overall.RData")

