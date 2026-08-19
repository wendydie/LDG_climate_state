# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 04_LDG_completeness_estimate.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(patchwork)
source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")
# Read the dataset----------------------------------------------------------
rich_df <- read.csv(sprintf("./results/LDG/%s_cell_%s_richness.csv", 
                            params$spacing, params$level))
time_bins <- readRDS("./data/time_bins.RDS")  # Load time bin information
lat_bins <- palaeoverse::lat_bins_area(n = 12) %>% arrange(min)

# Step 1: Data Processing --------------------------------------------------
# Only keep data where bin_midpoint < 486.85 (to match climate state data)
rich_df$stage <- time_bins$interval_name[match(rich_df$bin_midpoint, time_bins$mid_ma)]
rich_df <- rich_df %>%
  filter(bin_midpoint <= 486.8500) %>%
  mutate(
    completeness = ifelse(nT >= 5 & t <= 2 * nT, "Complete", "Incomplete"),
    abs_lat = abs(cell_lat),
    hemisphere = ifelse(cell_lat >= 0, "Northern", "Southern"),
    lat_band_mid = (floor(abs_lat / 30) * 30) + 15
  )
# Create 12 equal-area latitude bands for LDG slope calculation
rich_df <- rich_df %>%
  mutate(bin_index = findInterval(cell_lat, vec = c(lat_bins$min, Inf)),
         bin = lat_bins$bin[bin_index]) %>%
  left_join(lat_bins %>% select(bin, lat_bin_mid = mid), by = "bin") %>%
  mutate(abs_lat_bin_mid = abs(lat_bin_mid))  # Midpoints of 12 equal-area latitude bands

# -----------------------------------------------------------------------
# Determine "good" or "bad" using only Complete data
# Baseline QC: occurrence5_k1_tropical_temperate
# -----------------------------------------------------------------------

rich_df_complete <- rich_df %>%
  filter(completeness == "Complete") %>%
  mutate(
    lat_zone = case_when(
      abs_lat_bin_mid < 30 ~ "tropical",
      abs_lat_bin_mid < 60 ~ "temperate",
      abs_lat_bin_mid <= 90 ~ "polar",
      TRUE ~ NA_character_
    )
  )

adjacent_df <- has_adjacent_bins(rich_df_complete, lat_bins) %>%
  distinct(bin_midpoint, hemisphere, label) %>%
  transmute(
    bin_midpoint,
    hemisphere,
    has_adjacent_tt = label == "good"
  )

qc_label_df <- rich_df_complete %>%
  group_by(bin_midpoint, hemisphere) %>%
  summarise(
    has_tropical = any(lat_zone == "tropical", na.rm = TRUE),
    has_temperate = any(lat_zone == "temperate", na.rm = TRUE),
    n_valid_lat_bins = n_distinct(abs_lat_bin_mid),
    .groups = "drop"
  ) %>%
  left_join(
    adjacent_df,
    by = c("bin_midpoint", "hemisphere")
  ) %>%
  mutate(
    has_adjacent_tt = coalesce(has_adjacent_tt, FALSE),
    label = ifelse(
      has_tropical & has_temperate & has_adjacent_tt,
      "good",
      "bad"
    )
  )

# Merge hemisphere-level good/bad labels back into the full dataset
rich_df2 <- rich_df %>%
  left_join(
    qc_label_df %>%
      select(bin_midpoint, hemisphere, label),
    by = c("bin_midpoint", "hemisphere")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label)
  )

# Reorder completeness so "Complete" is at the bottom
rich_df2$completeness <- factor(rich_df2$completeness, levels = c("Incomplete", "Complete"))
# Prepare text annotation data (only for bad hemisphere)
bad_labels <- rich_df2 %>%
  select(bin_midpoint, cell_lat, hemisphere, label) %>%
  mutate(
    text_x = ifelse(hemisphere == "Northern", 45, -45),
    text_y = Inf  # Use Inf to place label at the top of each facet
  ) %>%
  filter(label == "bad") %>%
  distinct(bin_midpoint, hemisphere, text_x, text_y, label)

# Define hemisphere-specific background color
rich_df2 <- rich_df2 %>%
  mutate(bg_color = ifelse(label == "bad", "gray90", NA))  # Mark bad hemisphere

# Step 3: Plot Stacked Proportional Bar Chart ------------------------------
stacked_bar_plot <- ggplot(rich_df2, aes(x = cell_lat, fill = completeness)) +
  # Add a background rectangle for "bad" hemisphere bins
  geom_rect(data = rich_df2 %>% filter(label == "bad" & hemisphere == "Northern"), 
            aes(xmin = 0, xmax = 90, ymin = -Inf, ymax = Inf), 
            fill = "gray90", alpha = 0.3, inherit.aes = FALSE) +  
  geom_rect(data = rich_df2 %>% filter(label == "bad" & hemisphere == "Southern"), 
            aes(xmin = -90, xmax = 0, ymin = -Inf, ymax = Inf), 
            fill = "gray90", alpha = 0.3, inherit.aes = FALSE) +  
  # Add "Bad" text labels in bad hemispheres
  geom_text(data = bad_labels, aes(x = text_x, y = text_y, label = "Bad"), 
            color = "red", size = 3, fontface = "bold", vjust =1.5, inherit.aes = FALSE) +
  # Histogram
  geom_histogram(binwidth = 5, position = "stack", color = "black", linewidth = 0.3) +  
  # Define colors for Complete and Incomplete data
  scale_fill_manual(values = c("Incomplete" = "#F0E442","Complete" = "#009E73" )) +  
  # Equator & Baseline Lines
  geom_vline(xintercept = 0, color = "red", linetype = "solid", linewidth = 1) +  
  # Axis Labels
  labs(x = "Palaeolatitude", y = "Cell number", fill = "Completeness") +
  # Facet by 'stage' (bin_midpoint), sorted in descending order
  facet_wrap(~ reorder(bin_midpoint, -as.numeric(as.character(bin_midpoint))),
             labeller = as_labeller(function(x) 
               paste0(rich_df2$stage[match(x, rich_df2$bin_midpoint)])),
             nrow = 15, ncol = 6, scales = "free_y") +  
  # Keep x-axis strictly within [-90, 90] with ticks every 30 degrees
  scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) +  
  # Ensure Y-axis always has 4 labels (0, 1/3 max, 2/3 max, max)
  scale_y_continuous(
    breaks = function(y) {
      max_val <- ceiling(max(y, na.rm = TRUE) / 10) * 10  # Round up to nearest 10
      breaks <- seq(0, max_val, length.out = 5)  # Generate 5 evenly spaced breaks
      return(unique(round(breaks)))  # Ensure all breaks are integers
    }
  ) +
  # Improve theme appearance
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "bold", margin = margin(1,1,1,1)),  
    strip.placement = "inside",  
    panel.spacing = unit(0.01, "lines"),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    panel.grid = element_blank(),  
    axis.text.x = element_text(size = 8, color = "black"),  
    axis.title.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_text(size = 8, color = "black"),  
    axis.title.y = element_text(size = 12, color = "black", angle = 90),
    panel.background = element_rect(fill = "white"),
    # Adjust legend position and spacing
    legend.position = 'bottom',
    legend.margin = margin(-5, 0, 0, 0, unit = "pt"),  # Reduce top margin
    legend.spacing.y = unit(-3, "pt"),  # Reduce spacing between legend items
    legend.key.height = unit(5, "pt")  # Make legend more compact
  )
p_main <- stacked_bar_plot + theme(legend.position = "none")
g_legend <- cowplot::get_legend(
  stacked_bar_plot + 
    theme(
      legend.position = "bottom",
      legend.title = element_text(size=10),
      legend.text = element_text(size=8),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA), 
      legend.box.background = element_rect(fill = NA, color=NA)
    )
)
p_final <- p_main +
  inset_element(
    g_legend,
    left   = 0.63,
    bottom = 0.04,
    right  = 0.88,
    top    = 0.08,
    on_top = TRUE,
    clip = FALSE,
    align_to = "plot"
  )
print(p_final)

# Save and display plot
# print(stacked_bar_plot)
ggsave(sprintf("./figures/%s_km_%s_quota_completeness_histogram.jpg", 
               params$spacing, params$level), p_final, 
       width = 8, height = 8, dpi = 300)
