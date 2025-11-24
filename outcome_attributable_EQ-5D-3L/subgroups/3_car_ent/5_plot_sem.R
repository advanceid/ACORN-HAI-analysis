# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 extrafont,
                 ggplot2,
                 scales,
                 patchwork,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Data
factor_load <- readRDS("data/attfbis_factor_load_3_sub.RData")
S_coef <- readRDS("data/attfbis_S_coef_3_sub.RData")

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
  
  
  S_coef[[i]] <- S_coef[[i]] %>%
    mutate(sig_label = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  tit <- list("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")
  
  
  # Coef plot
  p <- ggplot(S_coef[[i]], aes(y = rhs, x = est.std)) +
    geom_col(width = 0.7, color = NA, fill = "steelblue") +  
    geom_errorbar(aes(xmin = boot.ci.lower, xmax = boot.ci.upper), 
                  width = 0.2, color = "black") +
    geom_text(aes(x = boot.ci.upper, label = paste0(pvalue, sig_label)), 
              hjust = -0.2, vjust = 0.5, size = 3, 
              color = "black", family = "Times New Roman") +
    theme_minimal() +
    scale_x_continuous(limits = c(-0.15, 0.5), breaks = seq(-0.15, 0.4, by = 0.1)) +
    labs(x = ifelse(i == 3, "Strength of standardized effect", ""), 
         y = NULL, title = tit[[i]]) +
    theme(
      text = element_text(family = "Times New Roman", size = 10, color = "black"),
      axis.text.y = element_text(size = 10, family = "Times New Roman", color = "black"),
      axis.text.x = if (i == 3) element_text(size = 10, family = "Times New Roman", 
                                             color = "black") else element_blank(),  
      axis.title.x = element_text(size = 10, family = "Times New Roman", color = "black"),
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


p <- wrap_plots(plots_combined, ncol = 1, heights = c(0.3, 0.5, 0.3))


# Save
CairoPDF("output/pdf/eq_3_sub.pdf", width = 9, height = 7)
print(p)
dev.off()