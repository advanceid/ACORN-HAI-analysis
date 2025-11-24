# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 extrafont,
                 ggplot2,
                 stringr,
                 scales,
                 patchwork,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
paths <- c("1_car_aci", "2_thir_ent", "3_car_ent")

factor_load <- list()
S_coef <- list()

for (i in seq_along(paths)) {
  factor_load[[i]] <- readRDS(file.path(paths[i], "data", paste0("attfbis_factor_load_", i, ".RData")))
  S_coef[[i]] <- readRDS(file.path(paths[i], "data", paste0("attfbis_S_coef_", i, ".RData")))
}

# Preprocess and plot
plots_combined <- list()

for (i in seq_along(S_coef)) {
  S_coef[[i]]$rhs <- factor(S_coef[[i]]$rhs, levels = rev(unique(S_coef[[i]]$rhs)))
  
  # ---- Order factor_load by est.std (largest at the top) ----
  ord_by_est <- factor_load[[i]] %>%
    dplyr::arrange(est.std) %>%   # ascending: small → large
    dplyr::pull(rhs)
  
  factor_load[[i]] <- factor_load[[i]] %>%
    dplyr::mutate(rhs = factor(rhs, levels = ord_by_est))
  
  S_coef[[i]]$pvalue <- as.character(S_coef[[i]]$pvalue)
  S_coef[[i]]$pvalue[S_coef[[i]]$pvalue == "0.05"] <- "<0.05"
  
  # special cases overwrite
  S_coef[[i]]$pvalue_num <- suppressWarnings(as.numeric(as.character(S_coef[[i]]$pvalue)))
  S_coef[[i]]$pvalue_num[S_coef[[i]]$pvalue == "<0.001"] <- 0.0009
  S_coef[[i]]$pvalue_num[S_coef[[i]]$pvalue == "<0.05"]  <- 0.049
  
  S_coef[[i]] <- S_coef[[i]] %>%
    mutate(
      sig_label = case_when(
        !is.na(pvalue_num) & pvalue_num < 0.001 ~ "***",
        !is.na(pvalue_num) & pvalue_num < 0.01  ~ "**",
        !is.na(pvalue_num) & pvalue_num < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  tit <- list(
    expression("Carbapenem-resistant" ~ italic("Acinetobacter") ~ "spp."),
    "Third-generation cephalosporin-resistant Enterobacterales",
    "Carbapenem-resistant Enterobacterales"
  )
  
  
  # Coef plot
  p <- ggplot(S_coef[[i]], aes(y = rhs, x = est.std)) +
    geom_col(width = 0.7, color = NA, fill = "steelblue") +  
    geom_errorbar(aes(xmin = boot.ci.lower, xmax = boot.ci.upper), 
                  width = 0.2, color = "black") +
    geom_text(aes(x = boot.ci.upper, label = paste0(pvalue, sig_label)), 
              hjust = -0.2, vjust = 0.5, size = 3, 
              color = "black", family = "Times New Roman") +
    theme_minimal() +
    scale_x_continuous(limits = c(-0.1, 0.42), breaks = seq(-0.1, 0.35, by = 0.1)) +
    labs(x = ifelse(i == 3, "Strength of standardized effect", ""), 
         y = NULL, title = tit[[i]]) +
    theme(
      text = element_text(family = "Times New Roman", size = 10, color = "black"),
      axis.text.y = element_text(size = 10, family = "Times New Roman", color = "black"),
      axis.text.x = if (i == 3) element_text(size = 10, family = "Times New Roman", 
                                             color = "black") else element_blank(),  
      axis.title.x = element_text(size = 10, family = "Times New Roman", color = "black",
                                  margin = margin(t = 10)),
      plot.title = element_text(size = 10, color = "black", 
                                family = "Times New Roman", hjust = 0.5),
      legend.title = element_text(size = 10, color = "black"), 
      legend.text = element_text(size = 10, color = "black", 
                                 margin = margin(l = 3, r = 6, unit = "pt")),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(t = -5, r = 0, b = 0, l = 0)
    )
  
  
  # Factor loading plot
  ann <- if (i == 1) {
    annotate("text", x = "EQ-5D-3L", y = length(factor_load[[i]]$rhs) + 1.2, 
             label = "EQ-5D-3L", size = 3.2, family = "Times New Roman")
  } else {
    NULL
  }
  
  p2 <- ggplot(factor_load[[i]], aes(y = rhs, x = "EQ-5D-3L", fill = est.std)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(rhs, " (", sprintf("%.3f", est.std), ")")),  
              color = "white", size = 3.2, fontface = "bold", family = "Times New Roman") +
    ann +
    expand_limits(y = length(factor_load[[i]]$rhs) + 2) +
    scale_fill_gradientn(colors = viridis_pal(option = "plasma", end = 0.9, direction = -1)(5),
                         name = "")  +
    theme_minimal() +
    labs(x = "", y = "", title = "") +  
    theme(
      text = element_text(family = "Times New Roman", size = 10, color = "black"),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  # Combine two plots
  plots_combined[[i]] <- p + p2 + plot_layout(ncol = 2, widths = c(1, 0.6))
}


p <- wrap_plots(plots_combined, ncol = 1)

# Save
CairoPDF("pdf/eq_overall.pdf", width = 9, height = 7)
print(p)
dev.off()
