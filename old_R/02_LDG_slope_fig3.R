# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_fig3.R
# Last updated: 2025-11-11
# -----------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)
# -----------------------------------------------------------------------
# Ensure consistent color mapping
rich_df <- rich_df %>%
  mutate(hemisphere_mod = factor(hemisphere_mod, levels = c("Northern", "Southern", "Poor quality")))
# Get unique colors present in the dataset
color_levels <- c("Northern", "Southern", "Poor quality")
color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Poor quality" = "#D3D3D3")
color_palette <- color_palette[color_levels]  # Keep only the necessary colors
# Define linetypes for different percentiles
shape_palette<- c("q50"=16,"q60"=17,"q75"=15,"q90"=18,"q95"=19)
quantile_palette <- c("q50" = "solid", "q60" = "dashed", "q75" = "dotdash", "q90" = "twodash", "q95" = "dotted")
# Create a lookup table between bin_midpoint and stage
lab_map <- rich_df %>%
  distinct(bin_midpoint, stage) %>%
  arrange(bin_midpoint)
stage_labeller <- setNames(lab_map$stage, lab_map$bin_midpoint)

# === 3) Split data into Northern and Southern Hemisphere subsets === #
# Also reshape percentile data (q50–q95) into long format for plotting
north_rich <- rich_df %>% filter(cell_lat > 0)
south_rich <- rich_df %>% filter(cell_lat < 0)

north_long <- north_rich %>%
  select(bin_midpoint, abs_lat_bin_mid, q50, q60, q75, q90, q95, hemisphere_mod) %>%
  pivot_longer(cols = c(q50, q60, q75, q90, q95),
               names_to = "quantile", values_to = "richness_value")

south_long <- south_rich %>%
  select(bin_midpoint, abs_lat_bin_mid, q50, q60, q75, q90, q95, hemisphere_mod) %>%
  pivot_longer(cols = c(q50, q60, q75, q90, q95),
               names_to = "quantile", values_to = "richness_value")

# Compute fitted Theil–Sen regression lines for each quantile
north_lines <- north_rich %>%
  select(bin_midpoint, stage, abs_lat_bin_mid, hemisphere, label,
         hemisphere_mod, qD_normalized, q50, q60, q75, q90, q95) %>%
  pivot_longer(cols = c(q50, q60, q75, q90, q95),
               names_to = "quantile", values_to = "qD_value") %>%
  left_join(LDG_slope[, c("bin_midpoint", "hemisphere", "quantile", "slope", "intercept")],
            by = c("bin_midpoint", "hemisphere", "quantile")) %>%
  mutate(fitted_values = intercept + slope * abs_lat_bin_mid) %>%
  filter(!is.na(fitted_values))

south_lines <- south_rich %>%
  select(bin_midpoint, stage, abs_lat_bin_mid, hemisphere, label,
         hemisphere_mod, qD_normalized, q50, q60, q75, q90, q95) %>%
  pivot_longer(cols = c(q50, q60, q75, q90, q95),
               names_to = "quantile", values_to = "qD_value") %>%
  left_join(LDG_slope[, c("bin_midpoint", "hemisphere", "quantile", "slope", "intercept")],
            by = c("bin_midpoint", "hemisphere", "quantile")) %>%
  mutate(fitted_values = intercept + slope * abs_lat_bin_mid) %>%
  filter(!is.na(fitted_values))

# === 4) Define a reusable plotting function === #
make_hemi_plot <- function(rich_side, long_side, lines_side, hemi_title) {
  ggplot(rich_side, aes(x = abs_lat, y = qD_normalized, color = hemisphere_mod)) +
    geom_point(alpha = 0.7, size = 1) +
    # Add percentile raw data points
    geom_point(data = long_side,
               aes(x = abs_lat_bin_mid, y = richness_value, color = hemisphere_mod, shape = quantile),
               size = 1, inherit.aes = FALSE) +
    # Add Theil–Sen regression lines for each percentile
    geom_line(data = lines_side,
              aes(x = abs_lat_bin_mid, y = fitted_values, linetype = quantile, color = hemisphere_mod),
              linewidth = 1, inherit.aes = FALSE) +
    scale_color_manual(name = "Line", values = color_palette, drop = TRUE) +
    scale_linetype_manual(name = "Percentile", values = quantile_palette) +
    scale_shape_manual(name = "Percentile", values = shape_palette) +
    guides(
      color    = guide_legend(title = "Line", override.aes = list(shape = NA, linetype = "solid")),
      linetype = guide_legend(title = "Percentile"),
      shape = "none"
    ) +
    scale_y_continuous(breaks = pretty_breaks(3),limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(
      ~ reorder(bin_midpoint, -as.numeric(as.character(bin_midpoint))),
      labeller = as_labeller(stage_labeller),
      ncol = 6
    ) +
    labs(
      title = hemi_title,
      x = "Absolute paleolatitude (°)",
      y = "Normalized generic richness"
    ) +
    theme_minimal() +
    theme(
      strip.text      = element_text(size = 8, face = "bold", margin = margin(1, 1, 1, 1)),
      strip.placement = "inside",
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.01, "lines"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid      = element_blank(),
      axis.text.x     = element_text(size = 8, color = "black"),
      axis.title.x    = element_text(size = 12, color = "black"),
      axis.text.y     = element_text(size = 8, color = "black"),
      axis.title.y    = element_text(size = 12, color = "black", angle = 90),
      legend.position = "bottom"
    )
}

# === 5) Generate separate plots for Northern and Southern Hemispheres === #
north_fig <- make_hemi_plot(north_rich, north_long, north_lines, "Northern Hemisphere")
south_fig <- make_hemi_plot(south_rich, south_long, south_lines, "Southern Hemisphere")

# === 6) Save both figures === #
north_path <- sprintf("./figures/LDG_slope_facet_%s_km_%s_quota_%s_equal_area_bins_NORTH.jpg",
                      params$spacing, params$level, rich_params$n_lat_bins)
south_path <- sprintf("./figures/LDG_slope_facet_%s_km_%s_quota_%s_equal_area_bins_SOUTH.jpg",
                      params$spacing, params$level, rich_params$n_lat_bins)

ggsave(north_path, north_fig, width = 8, height = 9, dpi = 300)
ggsave(south_path, south_fig, width = 8, height = 9, dpi = 300)

pal <- c("Northern" = "#0072B2",
         "Southern" = "#E69F00")

linetypes <- c("q50"="solid","q60"="dashed","q75"="dotdash","q90"="twodash","q95"="dotted")
shapes    <- c("q50"=16,"q60"=17,"q75"=15,"q90"=18,"q95"=19)
for (stg in unique(rich_df$stage)) {
  
  for (hemi in c("Northern","Southern")) {
    df_bin <- rich_df %>%
      filter(stage == stg, hemisphere == hemi)
    
    if (nrow(df_bin) == 0) next
    
    bin <- unique(df_bin$bin_midpoint)
    bin <- if (length(bin) > 0) bin[1] else NA
    
    olsl_data_bin <- ols_lines %>%
      filter(stage == stg, hemisphere == hemi)
    
    df_long <- df_bin %>%
      select(abs_lat_bin_mid, q50, q60, q75, q90, q95) %>%
      pivot_longer(c(q50, q60, q75, q90, q95),
                   names_to = "quantile", values_to = "richness_value")
    
    p <- ggplot(df_bin, aes(x = abs_lat, y = qD_normalized, color = hemisphere_mod)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_point(data = df_long,
                 aes(x = abs_lat_bin_mid, y = richness_value, shape = quantile),
                 size = 2, inherit.aes = FALSE, color = pal[hemi]) +
      geom_line(data = olsl_data_bin,
                aes(x = abs_lat_bin_mid, y = fitted_values, linetype = quantile),
                linewidth = 1, inherit.aes = FALSE, color = pal[hemi]) +
      
      # geom_smooth(method = "lm", se = FALSE, linetype = "dotted", linewidth = 1.1) +
      scale_color_manual(values = pal, guide = "none") +
      scale_linetype_manual(values = linetypes, name = "Percentile") +
      scale_shape_manual(values = shapes, name = "Percentile") +
      labs(
        title = sprintf("LDG Slope — %s (%s Ma) — %s", stg, bin, hemi),
        x = "Absolute paleolatitude (°)",
        y = "Normalized generic richness"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
    
    output_dir <- sprintf("./figures/LDG_slope_combined/%s km %squota %s equal_area latitude bins",
                          params$spacing, params$level, rich_params$n_lat_bins)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s_%s.jpg",
                         output_dir, bin, toupper(hemi))
    ggsave(file_name, p, width = 7, height = 5, dpi = 200)
  }
}
