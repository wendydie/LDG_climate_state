# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_sensitivity_test.R
# Last updated: 2025-10-15
# Author: Die (Wendy) Wen
# Email: geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)
# -----------------------------------------------------------------------
rich_df_test <- rich_df %>%
  group_by(bin_midpoint, hemisphere, abs_lat_bin_mid) %>%
  mutate(
    q50_raw = quantile(qD, 0.50, na.rm = TRUE),
    q60_raw = quantile(qD, 0.60, na.rm = TRUE),
    q75_raw = quantile(qD, 0.75, na.rm = TRUE),
    q90_raw = quantile(qD, 0.90, na.rm = TRUE),
    q95_raw = quantile(qD, 0.95, na.rm = TRUE),
    q50 = quantile(qD_normalized, 0.50, na.rm = TRUE),
    q60 = quantile(qD_normalized, 0.60, na.rm = TRUE),
    q75 = quantile(qD_normalized, 0.75, na.rm = TRUE),
    q90 = quantile(qD_normalized, 0.90, na.rm = TRUE),
    q95 = quantile(qD_normalized, 0.95, na.rm = TRUE),
  ) %>%
  ungroup()
# Define percentiles to analyze
percentiles <- c("q50", "q60", "q75", "q90", "q95")
color_palette <- c("N" = "#0072B2",
                   "S" = "#E69F00")
# Compute slopes separately for qD_normalized and qD for Northern and Southern Hemispheres
compute_slopes <- function(df, richness_column) {
  results <- list()
  
  for (perc in percentiles) {
    col_name <- if (richness_column == "qD") paste0(perc, "_raw") else perc
    slope_data <- df %>%
      group_by(bin_midpoint, hemisphere) %>%
      filter(n() > 2, length(unique(abs_lat_bin_mid)) > 1) %>%  # Ensure enough data points
      summarise(
        model = list(lm(reformulate("abs_lat_bin_mid", response = col_name), data = cur_data())),
        .groups = "drop"
      ) %>%
      rowwise() %>%
      mutate(
        slope = coef(model)[2],  # Extract Theil-Sen slope
        intercept = coef(model)[1],  # Extract Theil-Sen intercept
        quantile = perc,  # Store percentile for reference
        richness_type = richness_column  # Store richness type (qD or qD_normalized)
      ) %>%
      select(-model)
    
    results[[perc]] <- slope_data
  }
  
  return(bind_rows(results))
}
# -----------------------------------------------------------------------
# Compute slopes for both qD_normalized and qD
slope_qD_normalized <- compute_slopes(rich_df_test, "qD_normalized")
slope_qD <- compute_slopes(rich_df_test, "qD")

# Combine the percentile (quantile) and richness type (qD or qD_normalized)
LDG_slope_all <- bind_rows(slope_qD_normalized, slope_qD) %>%
  mutate(slope_type = paste0(quantile, ifelse(richness_type == "qD_normalized", "", "_raw")))

# Convert `slope_type` into wide format for correlation matrix
slope_data_wide <- LDG_slope_all %>%
  select(bin_midpoint, hemisphere, slope_type, slope) %>%
  pivot_wider(names_from = slope_type, values_from = slope)

# Drop bin_midpoint column for correlation calculation
slope_matrix <- slope_data_wide %>% 
  select(-bin_midpoint) %>%
  mutate(hemisphere = ifelse(hemisphere == "Northern", "N", "S"))

# Create the correlation scatterplot matrix
correlation_plot <- ggpairs(
  slope_matrix, 
  columns = 2:ncol(slope_matrix),  # Exclude "hemisphere" from numeric calculations
  aes(color = hemisphere, fill = hemisphere, alpha = 0.9),  # Color points by hemisphere
  upper = list(continuous = wrap("cor", size = 2)),  # Display correlation coefficient in upper panel
  lower = list(continuous = wrap("points", alpha = 0.7, size = 1)),  # Scatter plots in lower panel
  diag = list(continuous = wrap("densityDiag"))  # Density plots on diagonal
) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +  # Remove default ggplot2 background
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Keep black borders for subplots
    strip.background = element_blank(),  # Remove background (border) for facet labels
    strip.text = element_text(size = 8, face = "bold"),  # Reduce facet label text size
    axis.text = element_text(size = 6),  # Reduce axis text size
    axis.title = element_text(size = 8),  # Reduce axis title size
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Ensure no title border
    plot.background = element_blank()  # Ensure no background box for the whole plot
  )
print(correlation_plot)
cor_path <- sprintf("./figures/test/%skm %squota %slatitude band correlation plot.jpg", 
                    params$spacing, params$level, rich_params$n_lat_bins)
# Save the figure
ggsave(cor_path, correlation_plot, width = 8, height = 8, dpi = 300)
