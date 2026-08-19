# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 03_LDG_compared_in_climate_state.R
# Last updated: 2025-10-15
# Author: Die (Wendy) Wen
# Email: geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(grid)
library(patchwork)
library(cowplot)
library(deeptime)
library(ggnewscale)
source("./R/options.R")

# Wilcoxon signed-rank test analysis in different climate states---------
# Read dataset-----------------------------------------------------------
LDG_slope <- read.csv(sprintf("./results/%skm %squota %s equal-area latitude bins LDG slope.csv", 
                              params$spacing, params$level, rich_params$n_lat_bins))
climate_states <- read.csv("./data/climate_states.csv")
time_bins <- readRDS("./data/time_bins.RDS")
# Order climate states and climate colors
climate_levels <- c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse")
climate_colors <- c(
  "Coldhouse" = "#005344",
  "Coolhouse" = "#007d65",
  "Transitional" = "#c8c7c7",
  "Warmhouse" = "#b57a51",
  "Hothouse" = "#95484b"
)
#-------------------------------------------------------------------------
# Step 1: Merge final_results and climate_states by bin_midpoint----------
# Filter and assign colors
climate_states <- climate_states %>%
  filter(bottom <= 486.8500 & top >= 0) %>%
  mutate(climate_color = climate_colors[climate_state])  # Direct mapping
time_bins <- time_bins %>% filter(min_ma <= 486.8500 & max_ma >= 0)

slope_cli_df <- LDG_slope %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint))) %>%
  full_join(climate_states, by = c("bin_midpoint" = "mid"))

slope_cli_path <- sprintf("./results/%skm %squota %s equal-area latitude bins LDG slope and climate states.csv", 
                     params$spacing, params$level, rich_params$n_lat_bins)
write.csv(slope_cli_df, slope_cli_path, row.names = FALSE)

slope_cli_df_filter <- slope_cli_df[(slope_cli_df$climate_state != ''&
                                     slope_cli_df$label != 'bad' 
                                     &
                                     slope_cli_df$quantile == 'q75' # used the 75% to do the Wilcoxon signed-rank test
                                     ), ]

# Step 2: Group by climate_state and perform Wilcoxon signed-rank test-----
# Get all unique pairs of climate states
climate_pairs <- combn(climate_levels, 2, simplify = FALSE)

# Whether the LDG slope is normally distributed
test_results <- map_df(climate_pairs, function(pair) {
  group1 <- slope_cli_df_filter %>% filter(climate_state == pair[1]) %>% pull(slope)
  group2 <- slope_cli_df_filter %>% filter(climate_state == pair[2]) %>% pull(slope)
  
  # Perform Shapiro-Wilk normality test
  shapiro1_p <- shapiro.test(group1)$p.value
  shapiro2_p <- shapiro.test(group2)$p.value
  
  # Check if both groups are normally distributed
  normal <- shapiro1_p > 0.05 & shapiro2_p > 0.05
  test_type <- ifelse(normal, "t-test", "Wilcoxon")
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    normal = normal,
    test_type = test_type,
    shapiro_p1 = shapiro1_p,
    shapiro_p2 = shapiro2_p
  )
})
print(test_results)
test_path <- sprintf("./results/%skm %squota %s equal-area latitude bins selection of climate states.csv", 
                     params$spacing, params$level, rich_params$n_lat_bins)
write.csv(test_results, test_path, row.names = FALSE)

# Wilcoxon test results
wil_results <- map_df(climate_pairs, function(pair) {
  group1 <- slope_cli_df_filter %>% filter(climate_state == pair[1]) %>% pull(slope)
  group2 <- slope_cli_df_filter %>% filter(climate_state == pair[2]) %>% pull(slope)
  
  wilcox_test <- wilcox.test(group1, group2,alternative = "less") # Directional test: test if group1 < group2
  
  # Compute median difference safely
  median1 <- median(group1, na.rm = TRUE)
  median2 <- median(group2, na.rm = TRUE)
  median_diff <- ifelse(!is.na(median1) & !is.na(median2), median1 - median2, NA)
  
  # Compute additional statistics
  n1 <- sum(!is.na(group1))
  n2 <- sum(!is.na(group2))
  iqr1 <- IQR(group1, na.rm = TRUE) # IQR Quartile
  iqr2 <- IQR(group2, na.rm = TRUE)
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    n1 = n1,
    n2 = n2,
    median1 = median1,
    median2 = median2,
    p_value = wilcox_test$p.value,
    median_diff = median_diff,
    w_statistic = wilcox_test$statistic,
    iqr1 = iqr1,
    iqr2 = iqr2
  )
}) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adjusted)

wilr_path <- sprintf("./results/%skm %squota %s equal-area latitude bins wilcoxon test pairs of climate states.csv", 
                     params$spacing, params$level, rich_params$n_lat_bins)
write.csv(wil_results, wilr_path, row.names = FALSE)

#-------------------------------------------------------------------------
# Drawing the slope figures of LDG of time-bins in one figure-------------
# -- 1. Prepare slope data for plotting ---------------------
slope_data <- slope_cli_df %>%
  filter(quantile == 'q75') %>%
  mutate(
    slope_type = case_when(
      hemisphere == "Northern" ~ "Northern",
      hemisphere == "Southern" ~ "Southern",
      TRUE ~ NA_character_  # Handle unexpected cases
    ),
    slope_value = ifelse(label == "bad", NA, slope)
  )

# -- 2. Define axis ranges and theme -----------------------
data(periods)
data(epochs)
period_boundaries <- unique(c(periods$max_age, periods$min_age))
epoch_boundaries <- unique(c(epochs$max_age, epochs$min_age))
major_boundaries <- periods$max_age

x_max_val <- max(time_bins$max_ma) 
x_min_val <- min(time_bins$min_ma) 

y_min_val <- min(slope_data$slope_value, na.rm = TRUE) * 1.3
y_max_val <- max(slope_data$slope_value, na.rm = TRUE) *1.3
x_text_pos <- 15
y_text_pos <- y_min_val *0.94  # system text 

# Create a bar chart dataframe for climate state colors
climate_legend <- data.frame(
  climate_state = factor(c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse"),
                         levels = c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse")),
  max_ma = seq(x_max_val / 2 + 30,x_max_val / 2 - 30, length.out = 5),  # X-axis values
  min_ma = seq(x_max_val / 2 + 15,x_max_val / 2 - 45, length.out = 5),
  y = rep(y_max_val, 5),  # Ensure placement above the main plot
  climate_color = c("#005344", "#007d65", "#c8c7c7", "#b57a51", "#95484b")  # Color mapping
)

# Create a separate dataset for "Cooler Climate" and "Warmer Climate" labels
climate_labels <- data.frame(
  x = c(x_max_val / 2 + 35, x_max_val / 2 - 50),
  y = c(y_max_val * 0.8, y_max_val * 0.8),
  label = c("Cooler climate", "Warmer climate"),
  hjust = c(1,0)
)
climate_shade_layer <- list(
  geom_rect(
    data = climate_states,
    aes(xmin = bottom, xmax = top,
        ymin = y_min_val, ymax = y_max_val * 1.1,
        fill = I(climate_states$climate_color)),
    inherit.aes = FALSE, alpha = 0.6, colour = NA
  )
)
# -- 3. Draw the figure ----------------------------------
P_north_backcolor <- ggplot(filter(slope_data, slope_type == "Northern"), 
                  aes(x = bin_midpoint, y = slope_value, color = slope_type)) +
  climate_shade_layer +
  geom_vline(
    xintercept = major_boundaries,
    color = "black", linewidth = 0.4, alpha = 0.8
  ) +
  # Add climate state color bars
  geom_rect(
    data = climate_states, 
    aes(xmin = bottom, xmax = top, ymin = y_max_val, ymax = y_max_val * 1.1),
    fill = I(climate_states$climate_color), 
    color = "black", linewidth = 0.3, inherit.aes = FALSE
  ) +
  # annotate("rect",
  #          xmin = -Inf, xmax = Inf,
  #          ymin = -0.1, ymax = 0.1,
  #          fill = "lightblue", alpha = 0.3) +
  geom_hline(yintercept = 0, color="black", linewidth=0.4, linetype = "dashed") +
  # Add slope curves for different types
  geom_line(linewidth = 1, color = "#0072B2") +
  geom_point(
    aes(shape = abs(slope_value) < 0.1),
    size = 2, stroke = 0.55, color = "black", fill = "#0072B2"
  ) +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 21)) +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val*1.1), 
            color = "black", fill = NA, linewidth = 1) +
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val*1.1),
    breaks = seq(-3, 2, 1),
    expand = c(0, 0)
  )  +
  labs(
    x = NULL,
    y = 'Slope value',
    tag = "A",
    subtitle = "Northern Hemisphere"
  )+
  # annotate ("text", x = 540, y = y_max_val*0.85,label="A",
  #           size = 4, fontface = "bold")+
  # annotate ("text", x = x_text_pos, y = y_max_val*0.85,label="Northern Hemisphere",
  #           size = 4, fontface = "bold", hjust = 1)+
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = 'none',  # Remove redundant legend
    legend.title = element_text(size = 12, face = "bold"),  # Customize legend title
    legend.text = element_text(size = 10)  # Customize legend text
  )
P_south_backcolor <- ggplot(filter(slope_data, slope_type == "Southern"), 
                  aes(x = bin_midpoint, y = slope_value, color = slope_type)) +
  climate_shade_layer +
  geom_vline(
    xintercept = major_boundaries,
    color = "black", linewidth = 0.4, alpha = 0.8
  ) +
  # annotate("rect",
  #          xmin = -Inf, xmax = Inf,
  #          ymin = -0.1, ymax = 0.1,
  #          fill = "lightblue", alpha = 0.3) +
  geom_hline(yintercept = 0, color="black", linewidth=0.4, linetype = "dashed") +
  # Add slope curves for different types
  geom_line(linewidth = 1, color = "#E69F00") +
  geom_point(
    aes(shape = abs(slope_value) < 0.1),
    size = 2, stroke = 0.55, color = "black", fill = "#E69F00"
  ) +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 21)) +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val*1.1), 
            color = "black", fill = NA, linewidth = 1) +
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val*1.1),
    breaks = seq(-3, 2, 1),
    expand = c(0, 0)
  )  +
  labs(
    x = "Time (Ma)",
    y = "Slope value",
    tag = "B",
    subtitle = "Southern Hemisphere"
  ) +
  # annotate ("text", x = x_text_pos, y = y_max_val*0.85,label="Southern Hemisphere",
  #           size = 4, fontface = "bold", hjust = 1)+
  coord_geo(
    xlim = c(x_max_val, 0),
    pos = "bottom",
    dat = list("periods", "epochs"),
    height = unit(1.5, "lines"),
    expand = FALSE
  )+ 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = "none"
  )

climate_bar <- ggplot(climate_legend, aes(x = 0, y = 0, fill = climate_state)) +
  # Invisible points to trigger the legend (not plotted on canvas)
  geom_point(shape = 22, size = 5, alpha = 0) +
  scale_fill_manual(
    values = climate_colors,
    name   = "Climate state",
    drop   = FALSE
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, shape = 22, size = 5, colour = NA) # show color blocks in legend
    )
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 12),
    plot.margin     = margin(0, 10, 0, 10)
  )
# Extract only the legend as a grob
climate_legend_grob <- cowplot::get_legend(climate_bar)
# -- 4. Combine plots and print ---------------------------
slope_vTime_plot_backcolor <- ((P_north_backcolor / P_south_backcolor) / climate_legend_grob) +
  plot_layout(heights = c(10, 10, 1)) &
  theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.spacing = unit(0, "cm")
  )
# Display the plot
print(slope_vTime_plot_backcolor)
sT_path_jpg <- sprintf("./figures/jpg/background color %skm %squota %s equal-area latitude bins time series.jpg", 
                  params$spacing, params$level, rich_params$n_lat_bins)
ggsave(sT_path_jpg, slope_vTime_plot_backcolor, width = 8, height = 7, dpi = 900, type = "cairo")

sT_path_pdf <- sprintf("./figures/pdf/background color %skm %squota %s equal-area latitude bins time series.pdf", 
                   params$spacing, params$level, rich_params$n_lat_bins)
ggsave(sT_path_pdf, slope_vTime_plot_backcolor, width = 8, height = 7, dpi = 300)
#-------------------------------------------------------------------------
# Drawing the boxplot of slope of LDG in different climate states---------
# Filter LDG slope
slope_data_filtered <- slope_data %>%
  drop_na() %>%
  filter(climate_state %in% climate_levels) %>%
  mutate(climate_state = factor(climate_state, levels = climate_levels))

# Compute the total sample count for each climate state (combining all slope types)
sample_counts <- slope_data_filtered %>%
  select(climate_state, bin_midpoint) %>%
  group_by(climate_state) %>%  # Group only by climate state
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", count))  # Format label for display

state_n_labels <- sample_counts %>%
  mutate(
    state_with_n = paste0(climate_state, "\n(n=", count, ")"),
    climate_state = factor(climate_state, levels = climate_levels)
  ) %>%
  arrange(climate_state)

slope_data_filtered <- slope_data_filtered %>%
  left_join(state_n_labels[, c("climate_state", "state_with_n")], by = "climate_state") %>%
  mutate(state_with_n = factor(state_with_n, levels = state_n_labels$state_with_n))

slope_flag <- slope_data_filtered %>%
  group_by(state_with_n, slope_type) %>%
  mutate(
    q1 = quantile(slope_value, 0.25, na.rm = TRUE),
    q3 = quantile(slope_value, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr,
    is_outlier = slope_value < lower | slope_value > upper
  ) %>% 
  ungroup()

outlier_data <- slope_flag %>%
  mutate(slope_value = ifelse(is_outlier, slope_value, NA))

# Maximum x-axis position (for the rightmost arrow placement)
max_x <- max(as.numeric(factor(slope_data_filtered$climate_state))) + 1  

# Define y-axis boundaries
x_min_val <- min(as.numeric(factor(slope_data_filtered$climate_state))) - 0.5
x_max_val <- max_x - 0.5  
y_min_val <- min(slope_data_filtered$slope_value, na.rm = TRUE)
y_max_val <- max(slope_data_filtered$slope_value, na.rm = TRUE)

# **Boxplot layer**
boxplot <- ggplot(slope_flag, aes(x = state_with_n, y = slope_value, fill = slope_type)) +
  # Fill the "near zero" range from -0.1 to 0.1 with a light blue rectangle
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.1, ymax = 0.1,
           fill = "lightblue", alpha = 0.3) +
  # Boxplot layer
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8,linetype = "dashed") + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  # Jittered points layer - ensuring they are placed within the correct box
  # geom_jitter(aes(fill = slope_type),  # Ensure jitter uses correct fill colors
  #             # color = "gray",  # Use black border for points to differentiate them
  #             shape=21,
  #             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
  #             size = 1, alpha = 0.6) + 
  geom_jitter(
    data = subset(slope_flag, !is_outlier),
    shape = 21, size = 1, alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    show.legend = FALSE
  ) +
  geom_point(
    data = outlier_data,
    shape = 23, size = 1, stroke = 0.6, color = "black",
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  # Customizing colors for boxplot and jitter points
  scale_fill_manual(
    values = c("Northern" = "#0072B2",  # Blue for Northern slope
               "Southern" = "#E69F00"),  # Orange for Southern slope
    labels = c("Northern", "Southern")
  ) +
  
  # Labels and titles
  labs(
    x = "Climate state", y = "Slope value", fill = "Hemisphere"
  )  +
  # Minimal theme with some customizations
  theme_minimal() +                                       
  theme(axis.title = element_text(size = 14),                            
        axis.text = element_text(size = 12),
        axis.ticks = element_line(color = "black", linewidth=0.6),
        legend.title = element_text(size = 12),                         
        legend.text = element_text(size = 10),                         
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = c(0.02, 0.98),  # Move legend to top-left inside the plot
        legend.justification = c(0, 1),  # Align legend's top-left corner
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),  # Add a background
        legend.key = element_rect(fill = "white"),  # Keep legend keys clean
        plot.margin = margin(10, 30, 10, 10)
  ) +
  coord_cartesian(clip = 'off',xlim = c(1, 5), ylim = c(y_min_val,y_max_val)) +
  annotate("text", x = 5.81, y = (y_max_val + 0.1) / 2, label = "Non-modern-type", 
           vjust = 0.5, hjust = 0.5, size = 4.5, angle = 270) +
  annotate("text", x = 5.81, y = (y_min_val - 0.1) / 2, label = "Modern-type", 
           vjust = 0.5, hjust = 0.5, size = 4.5, angle = 270)

# Save and display
print(boxplot)

# Save high-resolution versions
bT_path_jpg <- sprintf("./figures/jpg/%skm %squota %s equal-area latitude bins boxplot.jpg", 
                       params$spacing, params$level, rich_params$n_lat_bins)
bT_path_pdf <- sprintf("./figures/pdf/%skm %squota %s equal-area latitude bins boxplot.pdf", 
                       params$spacing, params$level, rich_params$n_lat_bins)
ggsave(bT_path_jpg, boxplot, width = 7, height = 5, dpi = 300)
ggsave(bT_path_pdf, boxplot, width = 7, height = 5, dpi = 300)

