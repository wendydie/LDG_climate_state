# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_fig2.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(tidyr)
# -----------------------------------------------------------------------
# Convert q50, q60, q75, q90, q95 into long format for faceted plotting
df_long <- rich_df %>%
  select(bin_midpoint, abs_lat_bin_mid, q50, q60, q75, q90, q95, hemisphere_mod) %>%
  pivot_longer(cols = c(q50, q60, q75, q90, q95), 
               names_to = "quantile", values_to = "richness_value")
# Ensure consistent color mapping
rich_df <- rich_df %>%
  mutate(hemisphere_mod = factor(hemisphere_mod, levels = c("Northern", "Southern", "Poor quality")))
# Get unique colors present in the dataset
color_levels <- c("Northern", "Southern", "Poor quality")
color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Poor quality" = "#D3D3D3")
color_palette <- color_palette[color_levels]  # Keep only the necessary colors
# Define linetypes for different percentiles
quantile_palette <- c("q50" = "solid", "q60" = "dashed", "q75" = "dotdash", "q90" = "twodash", "q95" = "dotted")
# Create the faceted plot
combined_rich_fig <- ggplot(rich_df, aes(x = abs_lat, y = qD_normalized, color = hemisphere_mod)) +
  geom_point(alpha = 0.7, size = 1) +
  # Add raw richness percentile points
  geom_point(data = df_long, 
             aes(x = abs_lat_bin_mid, y = richness_value, color = hemisphere_mod, shape = quantile),
             size = 1, inherit.aes = FALSE) +
  # Add Theil-Sen fitted lines for each percentile
  geom_line(data = ols_lines, 
            aes(x = abs_lat_bin_mid, y = fitted_values, 
                linetype = quantile, color = hemisphere_mod), 
            linewidth = 1, inherit.aes = FALSE) +
  # Use viridis color palette (same for points and lines)
  scale_color_manual(name = "Line",
                     values = c("Poor quality" = "#D3D3D3",
                                "Northern" = "#0072B2",
                                "Southern" = "#E69F00")) +
  scale_linetype_manual(name = "Percentile Fit", values = quantile_palette) +
  scale_shape_manual(name = "Percentile", values = c("q50" = 16, "q60" = 17, "q75" = 15, "q90" = 18, "q95" = 19)) +
  guides(
    color = guide_legend(title = "Line",
                         override.aes = list(
      shape = NA,
      size  = 2,
      linetype = c("solid", "solid", "solid")
    )),
    linetype = guide_legend(title = "Percentile"),
    shape = "none"
  ) +
  facet_wrap(~ reorder(bin_midpoint, -as.numeric(as.character(bin_midpoint))),
             labeller = as_labeller(function(x) paste0(rich_df$stage[match(x, rich_df$bin_midpoint)])),
             ncol = 6) +
  # Facet by bin_midpoint with 8 columns
  scale_y_continuous(
    breaks = function(y) {
      max_val <- ceiling(max(y, na.rm = TRUE) / 10) * 10
      mid_val <- ceiling(max_val / 2 / 10) * 10
      return(c(0, mid_val, max_val))
    }
  ) +
  # Labels
  labs(
    x = "Absolute paleolatitude (°)",
    y = "Normalized generic richness"
  ) +
  theme_minimal() +
  # Reduce spacing between facet plots
  theme(
    strip.text = element_text(size = 8, face = "bold", margin = margin(1,1,1,1)),
    strip.placement = "inside",
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.01, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),  # Show x-axis labels at the bottom
    axis.title.x = element_text(size = 12, color = "black"),  # Show x-axis title
    axis.text.y = element_text(size = 8, color = "black"),  # Show y-axis labels on the left
    axis.title.y = element_text(size = 12, color = "black", angle = 90),  # Show y-axis title
    legend.position = "bottom"  # Remove redundant legend
  )
print(combined_rich_fig)
# Save the faceted figure
com_path  <- sprintf("./figures/LDG_slope_facet_%s_km_%s_quota_%s_equal_area_latitude_bins_combined_figure.jpg", 
                       params$spacing, params$level,rich_params$n_lat_bins)
ggsave(com_path, combined_rich_fig ,  width = 8, height = 9, dpi = 300)


for (stg in unique(rich_df$stage)) {

  df_bin <- rich_df %>%
    filter(stage == stg, !is.na(hemisphere_mod)) %>%
    mutate(hemisphere_mod = factor(hemisphere_mod, levels = c("Northern", "Southern", "Poor quality")))
  bin <- unique(df_bin$bin_midpoint)
  olsl_data_bin <- ols_lines %>%
    filter(stage == stg, !is.na(hemisphere_mod)) %>%
    mutate(hemisphere_mod = factor(hemisphere_mod, levels = c("Northern", "Southern", "Poor quality")))
  # Convert q50, q60, q75, q90, q95 into long format
  df_long <- df_bin %>%
    select(abs_lat_bin_mid, q50, q60, q75, q90, q95, hemisphere_mod) %>%
    pivot_longer(cols = c(q50, q60, q75, q90, q95),
                 names_to = "quantile", values_to = "richness_value")
  color_levels <- unique(df_bin$hemisphere_mod)
  color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Poor quality" = "#D3D3D3")
  color_palette <- color_palette[color_levels]  # Keep only the necessary colors
  # Define linetypes for different percentiles
  quantile_palette <- c("q50" = "solid", "q60" = "dashed", "q75" = "dotdash", "q90" = "twodash", "q95" = "dotted")

  p <- ggplot(df_bin, aes(x = abs_lat, y = qD_normalized, color = hemisphere_mod)) +
    geom_point(alpha = 0.7, size = 2) +
    # Add raw richness percentile points
    geom_point(data = df_long,
               aes(x = abs_lat_bin_mid, y = richness_value, color = hemisphere_mod, shape = quantile),
               size = 2, inherit.aes = FALSE) +
    # Add Theil-Sen fitted lines for each percentile
    geom_line(data = olsl_data_bin,
              aes(x = abs_lat_bin_mid, y = fitted_values,
                  linetype = quantile, color = hemisphere_mod),
              linewidth = 1, inherit.aes = FALSE) +
    # Add OLS fit for qD_normalized
    geom_smooth(data = df_bin, aes(x = abs_lat, y = qD_normalized, color = hemisphere_mod),
                method = "lm", se = FALSE, linetype = "dotted", linewidth = 1.2) +
    scale_color_manual(name = "Line", values = color_palette) +
    scale_linetype_manual(name = "Percentile fit", values = quantile_palette) +
    scale_shape_manual(name = "Percentiles", values = c("q50" = 16, "q60" = 17, "q75" = 15, "q90" = 18, "q95" = 19)) +
    guides(
      color = guide_legend(title = "Line", override.aes = list(shape = NA, linetype = "solid")), 
      # guide_legend(override.aes = list(color = unname(color_palette))),
      linetype = guide_legend(title = "Percentiles"),
      shape = "none",
    ) +
    labs(
      title = sprintf("LDG Slope for %s (%sMa)", stg, bin),
      x = "Absolute paleolatitude (°)",
      y = "Normalized generic richness"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  # Save the figure in a new output directory
  output_dir <- sprintf("./figures/LDG_slope_combined/%s km %squota %s equal_area latitude bins",
                        params$spacing, params$level, rich_params$n_lat_bins)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s_combined.jpg", output_dir, bin)
  ggsave(file_name, p, width = 7, height = 5, dpi = 200)
}
