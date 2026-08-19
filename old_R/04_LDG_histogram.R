# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 04_LDG_histogram.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
source("./R/options.R")
# -----------------------------------------------------------------------
# Read the dataset
rich_df <- read.csv(sprintf("./results/LDG/%s_cell_%s_richness.csv", 
                            params$spacing, params$level))
time_bins <- readRDS("./data/time_bins.RDS")  # Load time bin information

# Step 1: Data Processing --------------------------------------------------
# Only keep data where bin_midpoint < 486.85 (to match climate state data)
rich_df$stage <- time_bins$interval_name[match(rich_df$bin_midpoint, time_bins$mid_ma)]
rich_df <- rich_df %>%
  filter(bin_midpoint < 486.8500) %>%
  mutate(bin_midpoint = factor(bin_midpoint, levels = sort(unique(bin_midpoint), decreasing = TRUE)))

# Classify cells as "Incomplete" or "Complete"
rich_df <- rich_df %>%
  mutate(classification = ifelse(nT >= 5 & t <= 2 * nT, "Incomplete", "Complete"))

# Convert latitude to absolute values and classify hemispheres
rich_df <- rich_df %>%
  mutate(abs_lat = abs(cell_lat))

# Step 2: Prepare Data for Stacked Bar Plot ---------------------------------
# Aggregate counts by latitude and classification
rich_summary <- rich_df %>%
  group_by(cell_lat, classification) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cell_lat) %>%
  mutate(proportion = count / sum(count))  # Calculate proportion within each latitude bin

# Step 3: Plot Stacked Proportional Bar Chart ------------------------------
stacked_bar_plot <- ggplot(rich_summary, aes(x = cell_lat, y = proportion, fill = classification)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Stacked bar plot
  scale_fill_manual(values = c("Complete" = "#009E73", "Incomplete" = "#F0E442")) +  # Define colors
  scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, 30), expand = c(0, 0)) +  # Latitude axis ticks every 30 degrees
  labs(x = "Latitude (°)", y = "Proportion", fill = "Data Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Save and display plot
print(stacked_bar_plot)
ggsave(sprintf("./figures/%s_km_%s_quota_miltiton_rate_stacked_bar_plot.jpg", 
               params$spacing, params$level), stacked_bar_plot, width = 8, height = 6, dpi = 300)
