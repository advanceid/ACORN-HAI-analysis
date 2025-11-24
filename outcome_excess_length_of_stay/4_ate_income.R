# Clear workspace
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

# Read in the list of data frames for each cut-off
all_sub <- readRDS("data/los_lists.RData")

# Resistance pattern names
pattern_names <- c(
  "CRA vs CSA", "3GCRE vs 3GCSE", "CRE vs CSE",
  "CRP vs CSP", "VRE vs VSE", "MRSA vs MSSA"
)

# Cut-off labels
cutoffs <- c("t14", "t21", "t28", "t60", "t90")

# Stepwise AIC function
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

# Select covariates using t90
selected_covariates <- list()
df_list_t90 <- all_sub[["t90"]] 

for (i in seq_along(df_list_t90)) {
  df <- as.data.frame(df_list_t90[[i]])
  pat_name <- pattern_names[i]
  
  incomes <- levels(df$country_income) 
  
  for (inc in incomes) {
    df_sub <- df %>% filter(country_income == inc)
    
    candidate_vars <- c(
      "comorbidities_CCI",
      "severity_score_scale",
      "icu_hd_ap",
      "pathogen_combined_types"
    )
    
    selected_vars <- tryCatch({
      select_vars_logit(df_sub, "ast_ris", candidate_vars)
    }, error = function(e) character(0))
    
    # Always include age_new and sex
    selected_vars <- union(c("age_new", "sex"), selected_vars)
    selected_vars <- selected_vars[selected_vars != ""] 
    
    selected_covariates[[paste0(pat_name, "_", inc)]] <- selected_vars
    
    message(paste0("Selected (t90) for Pattern = ", pat_name, 
                   ", Income = ", inc, ": ", paste(selected_vars, collapse = ", ")))
  }
}

# Use fixed covariates for each cutoff
results_part1 <- list()

for (tc in cutoffs) {
  
  df_list <- all_sub[[tc]]
  t_int <- as.numeric(sub("t", "", tc))  
  
  for (i in seq_along(df_list)) {
    df <- as.data.frame(df_list[[i]])
    pat_name <- pattern_names[i]
    incomes <- levels(df$country_income)
    
    for (inc in incomes) {
      df_sub <- df %>% filter(country_income == inc)
      
      # Inclusion checks
      n_total <- nrow(df_sub)
      n_levels <- length(unique(df_sub$ast_ris))
      tab_event <- table(df_sub$event_new)
      n_evt1 <- ifelse("1" %in% names(tab_event), tab_event["1"], 0)
      
      n_analyzed <- n_total 
      
      if (n_total >= 100 && n_levels == 2 && n_evt1 >= 30) {
        
        # Use preselected variables
        selected_vars <- selected_covariates[[paste0(pat_name, "_", inc)]]
        if (is.null(selected_vars) || length(selected_vars) == 0) {
          # Use only a base formula if no valid variables were selected
          treat_formula <- as.formula("ast_ris ~ 1")
          selected_vars <- character(0) # For message consistency
        } else {
          treat_formula <- as.formula(paste("ast_ris ~", paste(selected_vars, collapse = " + ")))
        }
        
        # PS Calculation and 1%-99% Capping (Retain all data)
        ps_model <- tryCatch({
          glm(treat_formula, data = df_sub, family = binomial)
        }, error = function(e) {
          message(paste("PS Model Error:", e$message, "Skipping Capping for", pat_name, inc, tc))
          return(NULL)
        })
        
        if (!is.null(ps_model)) {
          
          ps_scores <- predict(ps_model, type = "response", na.action = na.pass)
          df_sub$ps <- ps_scores
          
          # Define Capping Thresholds
          trim_lower <- quantile(df_sub$ps, 0.025, na.rm = TRUE)
          trim_upper <- quantile(df_sub$ps, 0.975, na.rm = TRUE)
          
          # Apply Capping/Winsorizing to the PS score: replace extreme values
          df_sub$ps_capped <- pmin(pmax(df_sub$ps, trim_lower), trim_upper)
          
          # The analyzed N remains total N 
          n_analyzed <- n_total 
          n_capped <- sum(df_sub$ps != df_sub$ps_capped, na.rm = TRUE)
          
          message(paste0("Capping applied (", tc, ", ", pat_name, ", ", inc, 
                         "). PS values capped for: ", n_capped, 
                         " patients. Total N: ", n_analyzed, "."))
        } else {
          n_analyzed <- n_total
        }
        
        # cens_formula
        if (all(df_sub$event_new != 0)) {
          cens_formula <- ~1
        } else {
          cens_vars <- paste(c("ast_ris", selected_vars), collapse = " + ")
          cens_formula <- as.formula(paste("~", cens_vars))
        }
        
        
        # RMST model
        res1 <- tryCatch({
          fit <- resmeanATE(
            Event(time_new, event_new) ~ ast_ris,
            data = df_sub, # Full data used
            cause = 1,
            time = t_int,
            model = "exp",
            treat.model = treat_formula,
            augmentation = TRUE,
            cens.model = cens_formula,
            type = "II",
            prop.haz = FALSE,
            kaplan.meier = TRUE,
            se = TRUE
          )
          
          s <- summary(fit)$ateDR["treat:1-0", ]
          
          if (is.null(s) || all(is.na(s))) {
            message(paste("ATE estimate failed post-capping for:", pat_name, inc, tc))
            return(NULL)
          }
          
          data.frame(
            Cutoff = tc,
            Pattern = pat_name,
            Income = inc,
            Total_N = n_total,
            N_analyzed = n_analyzed, 
            Estimate = -round(as.numeric(s["Estimate"]), 0),
            CI_lower = -round(as.numeric(s["97.5%"]), 0),
            CI_upper = -round(as.numeric(s["2.5%"]), 0),
            P_value = as.numeric(s["P-value"]),
            stringsAsFactors = FALSE
          )
        }, error = function(e) {
          message(paste("Fatal Error in resmeanATE (Post-Capping) for:", pat_name, inc, tc, "-", e$message))
          return(NULL)
        })
        
        if (!is.null(res1)) {
          results_part1[[length(results_part1) + 1]] <- res1
        }
      }
    }
  }
}

# Bind and format
df_combined <- bind_rows(results_part1) %>%
  mutate(
    Group = paste0(Income, " (N = ", Total_N, ")"), 
    Estimate_CI = sprintf("%.0f [%.0f, %.0f]", Estimate, CI_lower, CI_upper),
    P_value = case_when(
      P_value < 0.001 ~ "<0.001",
      P_value < 0.01 ~ sprintf("%.3f", P_value),
      TRUE ~ sprintf("%.2f", P_value)
    )
  ) %>%
  select(Cutoff, Pattern, Group, Estimate_CI, P_value)

df_combined$group_by <- "World Bank income status"

# Save
write.csv(df_combined, "output/excess_los_by_income.csv", row.names = FALSE)
saveRDS(df_combined, "data/excess_los_by_income.RData")

