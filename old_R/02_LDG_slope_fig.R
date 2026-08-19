# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_fig.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(patchwork)
source("./R/options.R")
# -----------------------------------------------------------------------
rich_df <- rich_df %>%
  group_by(bin_midpoint, hemisphere, abs_lat_bin_mid) %>%
  mutate(
    q50 = quantile(qD_normalized, 0.50, na.rm = TRUE),
    q60 = quantile(qD_normalized, 0.60, na.rm = TRUE),
    q75 = quantile(qD_normalized, 0.75, na.rm = TRUE),
    q90 = quantile(qD_normalized, 0.90, na.rm = TRUE),
    q95 = quantile(qD_normalized, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()
rich_df <- rich_df %>%
  mutate(
    hemisphere_mod = ifelse(label == "bad", 
                            "Poor quality",
                            hemisphere)
  )
# Generate OLS fitted values for visualization, grouped by percentile
ols_lines <- rich_df %>%
  select(bin_midpoint, stage, abs_lat_bin_mid, hemisphere,label,hemisphere_mod,
         qD_normalized, q50, q60, q75, q90, q95) %>%  # Include percentiles
  pivot_longer(cols = c(q50, q60, q75, q90, q95), 
               names_to = "quantile", values_to = "qD_value") %>%  # Reshape to long format
  left_join(LDG_slope[,c("bin_midpoint", "hemisphere", "quantile", "slope", "intercept")], 
            by = c("bin_midpoint", "hemisphere", "quantile")) %>%  # Join with slopes based on bin, hemisphere, and quantile
  mutate(fitted_values = intercept + slope *  abs_lat_bin_mid) %>%  # Compute Theil-Sen fitted values
  filter(!is.na(fitted_values))  # Remove rows where fitted values could not be computed

# Step 4: Create the scatter plot with LDG slopes------------------------
LDG_s_plot <- ggplot(rich_df, aes(x = abs_lat, y = qD_normalized,
                                  color = hemisphere_mod,
                                  shape = hemisphere
)) +
  # Scatter plot points (Northern & Southern hemispheres get automatic colors)
  geom_point(alpha = 0.7, size = 1) +
  geom_point(data = filter(rich_df, label == "good"),
             aes(x = abs_lat_bin_mid, y = q75,shape = 'q75',
                 color = color),
             size = 2) +
  # Theil-Sen fitted line for Northern Hemisphere
  geom_line(data = filter(ols_lines, quantile=='q75'),
            aes(x = abs_lat_bin_mid, y = fitted_values, 
                linetype = hemisphere_mod, color = hemisphere_mod),
            linewidth = 1, inherit.aes = FALSE) +
  # Use viridis color palette (same for points and lines)
  scale_color_manual(name = "LDG slope",
                     values = c("Poor quality" = "#D3D3D3",
                                "Northern" = "#0072B2",
                                "Southern" = "#E69F00")) +
  scale_shape_manual(name = "Marker",
                     values = c("q75" = 4, "Northern" = 16, "Southern" = 17)) +
  scale_linetype_manual(name = "Line type", values = c("Poor quality" = "solid",
                                                       "Northern" = "solid",
                                                    "Southern" = "solid")) +
  # guides(
  #   shape = guide_legend(override.aes = list(size = c(2, 3, 3))),
  #   color = guide_legend(override.aes = list(color = c("#D3D3D3", "#0072B2", "#E69F00"), 
  #                                            shape = c(15, 16, 17))),
  #   linetype = "none"
  # ) +
  guides(
    shape = guide_legend(
      title = "Marker",direction = "horizontal",
      nrow = 1,
      override.aes = list(
        color = 'black',
        size  = c(2.5, 2, 2)
      )
    ),
    color = guide_legend(
      title = "Line",direction = "horizontal",
      nrow = 1,
      override.aes = list(
        shape = NA,
        size  = 2,
        linetype = c("solid", "solid", "solid")
      )
    ),
    linetype = "none"
  )+
  # Facet by bin_midpoint with 8 columns
  scale_y_continuous(
    breaks = function(y) {
      max_val <- ceiling(max(y, na.rm = TRUE) / 10) * 10
      mid_val <- ceiling(max_val / 2 / 10) * 10
      return(c(0, mid_val, max_val))
    }
  )  +
  facet_wrap(~ reorder(bin_midpoint, -as.numeric(as.character(bin_midpoint))),
             labeller = as_labeller(function(x) paste0(rich_df$stage[match(x, rich_df$bin_midpoint)])),
             ncol = 6) +
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
    panel.spacing = unit(0.01, "lines"),
    panel.spacing.x = unit(0.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),  # Show x-axis labels at the bottom
    axis.title.x = element_text(size = 12, color = "black"),  # Show x-axis title
    axis.text.y = element_text(size = 8, color = "black"),  # Show y-axis labels on the left
    axis.title.y = element_text(size = 12, color = "black", angle = 90),  # Show y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(-0.05,"pt")
  )
p_main <- LDG_s_plot + theme(legend.position = "none")
legend_obj <- cowplot::get_legend(
  LDG_s_plot +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(0,0,0,0),
      legend.margin = margin(0,0,0,0),
      legend.title = element_text(size=10),
      legend.text = element_text(size=8),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA), 
      legend.box.background = element_rect(fill = NA, color=NA)
    )
)
LDG_s_plot_final <- p_main +
  inset_element(
    legend_obj,
    left   = 0.63,
    bottom = -0.03,
    right  = 0.88,
    top    = 0.05,
    clip = FALSE,
    on_top = TRUE
  )
print(LDG_s_plot_final)
gg_path <- sprintf("./figures/%s km %squota %s equal-area latitude bins LDG slopes figure.jpg",
                   params$spacing, params$level, rich_params$n_lat_bins)
ggsave(gg_path, LDG_s_plot_final, width = 8, height = 9, dpi = 300)

for (stg in unique(rich_df$stage)) {
  df_bin <- subset(rich_df, stage == stg)
  bin <- unique(df_bin$bin_midpoint)
  for (perc in rich_params$percentiles) {
    olsl_data <- ols_lines %>% filter(stage == stg, quantile == perc)
    # Get unique color levels
    color_levels <- unique(df_bin$hemisphere_mod)
    # Define a fixed color palette, but only use the colors needed
    color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Poor quality" = "#D3D3D3")
    color_palette <- color_palette[color_levels]  # Keep only the required colors

    p <- ggplot(df_bin, aes(x = abs_lat, y = qD_normalized,
                            color = hemisphere_mod)) +  # Use new color column
      geom_point(alpha = 0.7, size = 2) +
      geom_point(data = filter(df_bin, label == "good"),
                 aes(x = abs_lat_bin_mid, y = get(perc),
                     color = hemisphere_mod), shape = 4, size = 2) +
      # Theil-Sen fitted lines (Merged for Northern & Southern Hemispheres)
      geom_line(data = olsl_data,
                aes(x = abs_lat_bin_mid, y = fitted_values, 
                    linetype = hemisphere, color = hemisphere_mod),
                linewidth = 1, inherit.aes = FALSE) +
      scale_color_manual(name = "LDG slope", values = color_palette) +  # Dynamically adjust colors
      scale_linetype_manual(name = "Legend", values = c("Northern" = "solid",
                                                        "Southern" = "solid")) +
      scale_x_continuous(
        limits = c(0, 90),
        expand = c(0, 0)
      )  +
      guides(
        color = guide_legend(override.aes = list(color = unname(color_palette))),
        linetype = "none"
      ) +
      labs(
        title = sprintf("LDG slope for %s (%s Ma)- percentile %s", stg, bin, perc),
        x = "Absolute paleolatitude (°)",
        y = sprintf("Generic richness (%s)", perc)
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
    # Save the figure
    output_dir <- sprintf("./figures/LDG slope per stage/%s km %squota %s equal_area latitude bins",
                          params$spacing, params$level, rich_params$n_lat_bins)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s_percentile_%s.jpg", output_dir, bin, perc)
    ggsave(file_name, p, width = 6, height = 5, dpi = 150)
  }
}

# # Bin latitude into intervals
# rich_df <- rich_df %>%
#   mutate(lat_band_mid = floor(cell_lat / rich_params$lat_band_width) * rich_params$lat_band_width)

# Compute richness percentiles across all bins
richness_percentiles <- rich_df %>%
  group_by(bin_midpoint, lat_bin_mid) %>%
  summarise(across(qD_normalized, list(
    q50 = ~quantile(., 0.50, na.rm = TRUE),
    q60 = ~quantile(., 0.60, na.rm = TRUE),
    q75 = ~quantile(., 0.75, na.rm = TRUE),
    q90 = ~quantile(., 0.90, na.rm = TRUE),
    q95 = ~quantile(., 0.95, na.rm = TRUE)
  ), .names = "{.fn}"), .groups = "drop") %>%
  pivot_longer(cols = starts_with("q"), names_to = "Percentile", values_to = "richness") %>%
  mutate(Percentile = factor(Percentile, levels = c("q50", "q60", "q75", "q90", "q95")))

# Use viridis color scheme (colorblind-friendly)
# percentile_colors <- setNames(viridis(5, option = "D"), c("q50", "q60", "q75", "q90", "q95"))
# percentile_colors <- setNames(gray.colors(5, start = 0.8, end=0.1),
#                               c("q50", "q60", "q75", "q90", "q95"))

# Plot species richness with faceting
LDG_fig <- ggplot() +
  geom_point(data = rich_df, aes(x = cell_lat, y = qD_normalized), size= 0.6, color = "black", alpha = 0.3) +  # Raw richness data
  geom_line(data = richness_percentiles, 
            aes(x = lat_bin_mid, y = richness, color = Percentile, group = Percentile),
            alpha = 0.6, linewidth = 0.8) +
  # scale_color_manual(name = "Percentile", values = percentile_colors) +
  labs(x = "Palaeolatitude",
       y = "Normalized generic richness") +
  # Facet by bin_midpoint with stage names
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
  )  +
  theme_minimal() +
  # Reduce spacing between facet plots
  theme(
    strip.text = element_text(size = 8, face = "bold", margin = margin(1,1,1,1)),
    strip.placement = "inside",
    panel.spacing.y = unit(0.01, "lines"),
    panel.spacing.x = unit(0.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),  # Show x-axis labels at the bottom
    axis.title.x = element_text(size = 12, color = "black"),  # Show x-axis title
    axis.text.y = element_text(size = 8, color = "black"),  # Show y-axis labels on the left
    axis.title.y = element_text(size = 12, color = "black", angle = 90),  # Show y-axis title
    legend.position = "bottom"
  )
p_fig_main <- LDG_fig + theme(legend.position = "none")
fig_legend_obj <- cowplot::get_legend(
  LDG_fig +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(0,0,0,0),
      legend.margin = margin(0,0,0,0),
      legend.title = element_text(size=10),
      legend.text = element_text(size=8),
      legend.spacing.x = unit(-0.05,"pt"),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA), 
      legend.box.background = element_rect(fill = NA, color=NA)
    )
)
LDG_fig_final <- p_fig_main +
  inset_element(
    fig_legend_obj,
    left   = 0.63,
    bottom = 0.01,
    right  = 0.88,
    top    = 0.03,
    clip = FALSE,
    on_top = TRUE
  )
print(LDG_fig_final)
# Save the faceted figure
LDG_f_path <- sprintf("./figures/LDG_per_stage_facet_%s_km_%s_quota %s_euqal_area_latitude_bins.jpg",
                       params$spacing, params$level,rich_params$n_lat_bins)
ggsave(LDG_f_path, LDG_fig_final, width = 8, height = 9, dpi = 300)


for (stg in unique(rich_df$stage)) {
  
  df_bin <- subset(rich_df, stage == stg)
  bin <- unique(df_bin$bin_midpoint)
  # Compute richness percentiles
  richness_percentile <- df_bin %>%
    group_by(lat_bin_mid) %>%
    summarise(across(qD_normalized, list(
      q50 = ~quantile(., 0.50, na.rm = TRUE),
      q60 = ~quantile(., 0.60, na.rm = TRUE),
      q75 = ~quantile(., 0.75, na.rm = TRUE),
      q90 = ~quantile(., 0.90, na.rm = TRUE),
      q95 = ~quantile(., 0.95, na.rm = TRUE)),
                     .names = "{.fn}")) %>%
    pivot_longer(cols = starts_with("q"), names_to = "Percentile", values_to = "richness") %>%
    ungroup() %>%
    mutate(Percentile = factor(Percentile, levels = c("q50", "q60", "q75", "q90", "q95")))  # Ensure correct legend order
  # Plot species richness
  p <- ggplot() +
    geom_point(data = df_bin, aes(x = cell_lat, y = qD_normalized), color = "black", alpha = 0.5) +  # Raw richness data
    geom_line(data = richness_percentile, aes(x = lat_bin_mid, y = richness, color = Percentile, group = Percentile), size = 1) +
    # scale_color_manual(name = "Percentile", values = percentile_colors) +
    scale_x_continuous(
      limits = c(-90, 90),
      expand = c(0, 0)
    )  +
    labs(title = sprintf("Genus Richness in %s (%s Ma)", stg, bin),
         x = "Palaeolatitude",
         y = "Normalized richness") +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  # Save the figure
  output_dir <- sprintf("./figures/LDG per stage/%s km %squota %s equal_area latitude bins",
                        params$spacing, params$level, rich_params$n_lat_bins)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s.jpg", output_dir, bin)
  ggsave(file_name, p, width = 6, height = 5, dpi = 150)
}

