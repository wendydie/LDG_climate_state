# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File: 02_03_LDG_slope_per_cell_allcells_climate_state.R
# Purpose:
#   Estimate LDG slopes using all individual grid-cell richness estimates,
#   compare slopes among climate states, and draw time-series and boxplot.
# -----------------------------------------------------------------------

# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(palaeoverse)
  library(deeptime)
  library(patchwork)
  library(cowplot)
  library(grid)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

# -----------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------

occurrence_min <- 5

climate_levels <- c(
  "Coldhouse",
  "Coolhouse",
  "Transitional",
  "Warmhouse",
  "Hothouse"
)

# climate_colors <- c(
#   "Coldhouse"    = "#005344",
#   "Coolhouse"    = "#007d65",
#   "Transitional" = "#c8c7c7",
#   "Warmhouse"    = "#b57a51",
#   "Hothouse"     = "#95484b"
# )
climate_colors <- c(
  "Coldhouse"    = "#B2CCC7",
  "Coolhouse"    = "#B3D8D0",
  "Transitional" = "#E4E4E4",
  "Warmhouse"    = "#DECABC",
  "Hothouse"     = "#D6B6B7"
)
hemi_cols <- c(
  "Northern" = "#0072B2",
  "Southern" = "#E69F00"
)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
dir.create("./figures/jpg", recursive = TRUE, showWarnings = FALSE)
dir.create("./figures/pdf", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 1. Helper functions
# -----------------------------------------------------------------------

fit_ols_safe <- function(df,
                         y_col = "qD_normalized",
                         x_col = "abs_lat") {
  
  df <- df %>%
    filter(
      !is.na(.data[[y_col]]),
      !is.na(.data[[x_col]])
    )
  
  if (nrow(df) < 3 || n_distinct(df[[x_col]]) < 2) {
    return(tibble(
      slope = NA_real_,
      intercept = NA_real_,
      slope_se = NA_real_,
      slope_lower_95 = NA_real_,
      slope_upper_95 = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_cells_model = nrow(df)
    ))
  }
  
  mod <- tryCatch(
    lm(as.formula(sprintf("%s ~ %s", y_col, x_col)), data = df),
    error = function(e) NULL
  )
  
  if (is.null(mod)) {
    return(tibble(
      slope = NA_real_,
      intercept = NA_real_,
      slope_se = NA_real_,
      slope_lower_95 = NA_real_,
      slope_upper_95 = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_cells_model = nrow(df)
    ))
  }
  
  sm <- summary(mod)
  cf <- sm$coefficients
  
  if (!(x_col %in% rownames(cf))) {
    return(tibble(
      slope = NA_real_,
      intercept = NA_real_,
      slope_se = NA_real_,
      slope_lower_95 = NA_real_,
      slope_upper_95 = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_cells_model = nrow(df)
    ))
  }
  
  slope <- as.numeric(cf[x_col, "Estimate"])
  slope_se <- as.numeric(cf[x_col, "Std. Error"])
  
  tibble(
    slope = slope,
    intercept = as.numeric(cf["(Intercept)", "Estimate"]),
    slope_se = slope_se,
    slope_lower_95 = slope - 1.96 * slope_se,
    slope_upper_95 = slope + 1.96 * slope_se,
    p_value = as.numeric(cf[x_col, "Pr(>|t|)"]),
    r_squared = as.numeric(sm$r.squared),
    n_cells_model = nrow(df)
  )
}

safe_shapiro <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  
  if (length(x) < 3 || length(unique(x)) < 3) {
    return(NA_real_)
  }
  
  shapiro.test(x)$p.value
}

safe_wilcox_less <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  if (length(x) < 1 || length(y) < 1) {
    return(list(p_value = NA_real_, statistic = NA_real_))
  }
  
  wt <- wilcox.test(x, y, alternative = "less", exact = FALSE)
  
  list(
    p_value = wt$p.value,
    statistic = as.numeric(wt$statistic)
  )
}

# -----------------------------------------------------------------------
# 2. Read and prepare data
# -----------------------------------------------------------------------

rich_df <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,
  params$level
))

time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(n = rich_params$n_lat_bins) %>%
  arrange(min) %>%
  mutate(
    abs_lat_bin_mid = round(abs(mid), 3),
    lat_zone = case_when(
      abs_lat_bin_mid < 30 ~ "Low",
      abs_lat_bin_mid < 60 ~ "Middle",
      abs_lat_bin_mid <= 90 ~ "High",
      TRUE ~ NA_character_
    )
  )

lat_zone_lookup <- lat_bins %>%
  select(bin, abs_lat_bin_mid, lat_zone)

rich_df <- rich_df %>%
  filter(
    nT >= occurrence_min,
    t <= 2 * nT
  ) %>%
  mutate(
    stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)]
  ) %>%
  filter(bin_midpoint <= 486.8500) %>%
  mutate(
    bin_index = findInterval(cell_lat, vec = c(lat_bins$min, Inf)),
    bin = lat_bins$bin[bin_index],
    abs_lat = abs(cell_lat),
    hemisphere = case_when(
      cell_lat >= 0 ~ "Northern",
      cell_lat < 0  ~ "Southern",
      TRUE ~ NA_character_
    ),
    lat_band_mid = floor(abs_lat / 30) * 30 + 15
  ) %>%
  left_join(lat_zone_lookup, by = "bin") %>%
  filter(
    !is.na(abs_lat),
    !is.na(hemisphere),
    !is.na(abs_lat_bin_mid),
    !is.na(lat_zone)
  ) %>%
  group_by(bin_midpoint) %>%
  mutate(
    qD_normalized = qD * 100 / max(qD, na.rm = TRUE)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------
# 3. Original QC labels unchanged
# -----------------------------------------------------------------------

if (rich_params$n_lat_bins == 6) {
  
  south_bin_no <- lat_bins$bin[lat_bins$mid >= -40 & lat_bins$mid <= 0]
  north_bin_no <- lat_bins$bin[lat_bins$mid <= 40 & lat_bins$mid >= 0]
  
  rich_df <- rich_df %>%
    group_by(bin_midpoint, hemisphere) %>%
    mutate(
      label = ifelse(
        (hemisphere == "Southern" & all(south_bin_no %in% bin)) |
          (hemisphere == "Northern" & all(north_bin_no %in% bin)),
        "good",
        "bad"
      )
    ) %>%
    ungroup()
  
} else {
  
  rich_df <- has_adjacent_bins(rich_df, lat_bins)
  
}

rich_df <- rich_df %>%
  mutate(color = ifelse(label == "bad", "Bad hemisphere", hemisphere))

# -----------------------------------------------------------------------
# 4. Calculate all-cell per-cell OLS slopes
# -----------------------------------------------------------------------

LDG_slope <- rich_df %>%
  group_by(bin_midpoint, stage, hemisphere) %>%
  group_modify(~ fit_ols_safe(
    .x,
    y_col = "qD_normalized",
    x_col = "abs_lat"
  )) %>%
  ungroup() %>%
  mutate(
    method = "grid_cell_all_cells_OLS"
  )

label_df <- rich_df %>%
  distinct(bin_midpoint, stage, hemisphere, label, color)

sampling_summary <- rich_df %>%
  group_by(bin_midpoint, stage, hemisphere) %>%
  summarise(
    n_cells = n(),
    n_valid_lat_bins = n_distinct(abs_lat_bin_mid),
    min_abs_lat = min(abs_lat, na.rm = TRUE),
    max_abs_lat = max(abs_lat, na.rm = TRUE),
    lat_range = max_abs_lat - min_abs_lat,
    sum_nT = sum(nT, na.rm = TRUE),
    mean_nT = mean(nT, na.rm = TRUE),
    median_nT = median(nT, na.rm = TRUE),
    mean_t_over_nT = mean(t / nT, na.rm = TRUE),
    max_t_over_nT = max(t / nT, na.rm = TRUE),
    prop_high_extrapolation = mean(t / nT > 1.5, na.rm = TRUE),
    .groups = "drop"
  )

LDG_slope <- LDG_slope %>%
  left_join(label_df, by = c("bin_midpoint", "stage", "hemisphere")) %>%
  left_join(sampling_summary, by = c("bin_midpoint", "stage", "hemisphere")) %>%
  mutate(
    slope_direction = case_when(
      slope < 0 ~ "normal_LDG",
      slope > 0 ~ "reverse_LDG",
      slope == 0 ~ "flat_LDG",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(hemisphere, bin_midpoint)

slope_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell all-cells OLS.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

write.csv(LDG_slope, slope_path, row.names = FALSE)

# -----------------------------------------------------------------------
# 5. Merge with climate states
# -----------------------------------------------------------------------

climate_states <- read.csv("./data/climate_states.csv") %>%
  filter(bottom <= 486.8500 & top >= 0) %>%
  mutate(
    climate_state = factor(climate_state, levels = climate_levels),
    climate_color = climate_colors[as.character(climate_state)]
  )

time_bins <- time_bins %>%
  filter(min_ma <= 486.8500 & max_ma >= 0)

slope_cli_df <- LDG_slope %>%
  mutate(
    bin_midpoint = as.numeric(as.character(bin_midpoint)),
    slope = as.numeric(as.character(slope))
  ) %>%
  full_join(climate_states, by = c("bin_midpoint" = "mid"))

slope_cli_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell all-cells OLS and climate states.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

write.csv(slope_cli_df, slope_cli_path, row.names = FALSE)

slope_cli_df_filter <- slope_cli_df %>%
  mutate(
    slope = as.numeric(as.character(slope)),
    climate_state = as.character(climate_state)
  ) %>%
  filter(
    climate_state != "",
    climate_state %in% climate_levels,
    label != "bad",
    !is.na(slope)
  )

# -----------------------------------------------------------------------
# 6. Climate-state Wilcoxon tests
# -----------------------------------------------------------------------

climate_pairs <- combn(climate_levels, 2, simplify = FALSE)

test_results <- map_df(climate_pairs, function(pair) {
  
  g1 <- slope_cli_df_filter %>%
    filter(climate_state == pair[1]) %>%
    pull(slope) %>%
    as.numeric()
  
  g2 <- slope_cli_df_filter %>%
    filter(climate_state == pair[2]) %>%
    pull(slope) %>%
    as.numeric()
  
  shapiro1_p <- safe_shapiro(g1)
  shapiro2_p <- safe_shapiro(g2)
  
  normal <- !is.na(shapiro1_p) &
    !is.na(shapiro2_p) &
    shapiro1_p > 0.05 &
    shapiro2_p > 0.05
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    normal = normal,
    test_type = ifelse(normal, "t-test", "Wilcoxon"),
    shapiro_p1 = shapiro1_p,
    shapiro_p2 = shapiro2_p
  )
})

wil_results <- map_df(climate_pairs, function(pair) {
  
  g1 <- slope_cli_df_filter %>%
    filter(climate_state == pair[1]) %>%
    pull(slope) %>%
    as.numeric()
  
  g2 <- slope_cli_df_filter %>%
    filter(climate_state == pair[2]) %>%
    pull(slope) %>%
    as.numeric()
  
  g1 <- g1[!is.na(g1)]
  g2 <- g2[!is.na(g2)]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  if (n1 < 1 || n2 < 1) {
    return(tibble(
      group1 = pair[1],
      group2 = pair[2],
      n1 = n1,
      n2 = n2,
      median1 = NA_real_,
      median2 = NA_real_,
      p_value = NA_real_,
      median_diff = NA_real_,
      w_statistic = NA_real_,
      iqr1 = NA_real_,
      iqr2 = NA_real_
    ))
  }
  
  wt <- safe_wilcox_less(g1, g2)
  
  med1 <- median(g1, na.rm = TRUE)
  med2 <- median(g2, na.rm = TRUE)
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    n1 = n1,
    n2 = n2,
    median1 = med1,
    median2 = med2,
    p_value = wt$p_value,
    median_diff = med1 - med2,
    w_statistic = wt$statistic,
    iqr1 = IQR(g1, na.rm = TRUE),
    iqr2 = IQR(g2, na.rm = TRUE)
  )
}) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adjusted)

test_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins per-cell all-cells OLS normality tests.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

wil_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins per-cell all-cells OLS wilcoxon tests.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

write.csv(test_results, test_path, row.names = FALSE)
write.csv(wil_results, wil_path, row.names = FALSE)

print(wil_results)

# -----------------------------------------------------------------------
# 7. Time-series figure
# -----------------------------------------------------------------------

data(periods)
data(epochs)

major_boundaries <- periods$max_age

slope_data <- slope_cli_df %>%
  mutate(
    slope = as.numeric(as.character(slope)),
    slope_type = case_when(
      hemisphere == "Northern" ~ "Northern",
      hemisphere == "Southern" ~ "Southern",
      TRUE ~ NA_character_
    ),
    uncertain = case_when(
      label == "bad" ~ NA,
      is.na(slope_lower_95) | is.na(slope_upper_95) ~ NA,
      slope_lower_95 <= 0 & slope_upper_95 >= 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    slope_value = ifelse(label == "bad", NA_real_, slope),
    slope_lower_plot = ifelse(label == "bad", NA_real_, slope_lower_95),
    slope_upper_plot = ifelse(label == "bad", NA_real_, slope_upper_95)
  )

x_max_val <- max(time_bins$max_ma, na.rm = TRUE)
x_min_val <- min(time_bins$min_ma, na.rm = TRUE)

major_boundaries_plot <- major_boundaries[
  major_boundaries >= 0 & major_boundaries <= x_max_val
]

y_min_val <- min(
  c(slope_data$slope_value, slope_data$slope_lower_plot),
  na.rm = TRUE
) * 1.01

y_max_val <- max(
  c(slope_data$slope_value, slope_data$slope_upper_plot),
  na.rm = TRUE
) * 1.1

if (
  y_min_val == y_max_val ||
  !is.finite(y_min_val) ||
  !is.finite(y_max_val)
) {
  y_min_val <- -1
  y_max_val <- 1
}

y_breaks <- pretty(c(y_min_val, y_max_val), n = 5)

climate_shade_layer <- list(
  geom_rect(
    data = climate_states,
    aes(
      xmin = bottom,
      xmax = top,
      ymin = y_min_val,
      ymax = y_max_val * 1.1,
      fill = I(climate_color)
    ),
    inherit.aes = FALSE,
    colour = NA
  )
)

plot_slope_ts <- function(
  df,
  hemi,
  col,
  tag_lab,
  show_top_bar = FALSE,
  show_x_tick = FALSE,
  bottom_space = 0
) {
  p <- ggplot(
    filter(df, slope_type == hemi),
    aes(x = bin_midpoint, y = slope_value)
  ) +
    climate_shade_layer +
    geom_vline(
      xintercept = major_boundaries_plot,
      color = "black",
      linewidth = 0.35,
      alpha = 0.7
    )

  if (show_top_bar) {
    p <- p +
      geom_rect(
        data = climate_states,
        aes(
          xmin = bottom,
          xmax = top,
          ymin = y_max_val,
          ymax = y_max_val * 1.1
        ),
        fill = I(climate_states$climate_color),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
  }

  p <- p +
    geom_hline(
      yintercept = 0,
      color = "black",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    geom_linerange(
      aes(
        ymin = slope_lower_plot,
        ymax = slope_upper_plot
      ),
      colour = unname(col),
      alpha = 0.8,
      linewidth = 0.75,
      na.rm = TRUE
    ) +
    geom_line(
      linewidth = 1,
      color = unname(col),
      na.rm = TRUE
    ) +
    geom_point(
      shape=21,
      size = 2,
      stroke = 0.5,
      color = "black",
      fill = unname(col),
      na.rm = TRUE
    ) +
    # scale_shape_manual(
    #   values = c(`TRUE` = 1, `FALSE` = 21),
    #   na.translate = FALSE
    # ) +
    annotate(
      "rect",
      xmin = x_min_val,
      xmax = x_max_val,
      ymin = y_min_val,
      ymax = y_max_val * 1.1,
      fill = NA,
      color = "black",
      linewidth = 0.8
    ) +
    scale_x_reverse(
      name = NULL,
      limits = c(x_max_val, 0),
      breaks = seq(500, 0, -50),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(y_min_val, y_max_val * 1.1),
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Slope value",
      tag = tag_lab,
      subtitle = paste(hemi, "Hemisphere")
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      plot.subtitle = element_text(size = 12, face = "bold"),
      plot.tag = element_text(size = 13, face = "bold"),
      plot.tag.position = c(0.01, 0.98),
      legend.position = "none",
      plot.margin = margin(0, 0, bottom_space, 0)
    )

  if (show_x_tick) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5)
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  p
}

P_north_backcolor <- plot_slope_ts(
  slope_data,
  "Northern",
  hemi_cols["Northern"],
  "(a)",
  show_top_bar = TRUE,
  show_x_tick = TRUE,
  bottom_space = 8
)

P_south_backcolor <- plot_slope_ts(
  slope_data,
  "Southern",
  hemi_cols["Southern"],
  "(b)",
  show_top_bar = FALSE,
  show_x_tick = FALSE,
  bottom_space = 0
)

P_geo <- ggplot(
  data.frame(x = c(x_max_val, 0), y = c(0, 0)),
  aes(x = x, y = y)
) +
  geom_blank() +
  scale_x_reverse(
    name = "Time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 0.01), expand = c(0, 0)) +
  coord_geo(
    xlim = c(x_max_val, 0),
    pos = "bottom",
    dat = list("periods", "epochs"),
    height = unit(1.45, "lines"),
    expand = FALSE
  ) +
  labs(x = "Time (Ma)", y = NULL, tag = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 14, margin = margin(t = 2)),
    axis.text.x = element_text(size = 14),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.tag = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 6, 0)
  )

climate_bar <- ggplot(
  data.frame(climate_state = factor(climate_levels, levels = climate_levels)),
  aes(x = climate_state, y = 1, fill = climate_state)
) +
  geom_tile(width = 0.8, height = 0.8) +
  scale_fill_manual(
    values = climate_colors,
    name = "Climate state",
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 1, colour = NA)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.margin = margin(4, 0, 4, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  )

climate_legend_grob <- cowplot::get_legend(climate_bar)

geo_climate_spacer <- cowplot::ggdraw()

aligned_plots <- cowplot::align_plots(
  P_north_backcolor,
  P_south_backcolor,
  P_geo,
  align = "v",
  axis = "lr"
)

slope_vTime_plot <- cowplot::plot_grid(
  aligned_plots[[1]],
  aligned_plots[[2]],
  aligned_plots[[3]],
  geo_climate_spacer,
  climate_legend_grob,
  ncol = 1,
  rel_heights = c(10, 10, 1.65, 0.45, 1.15)
) +
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )

print(slope_vTime_plot)

ts_jpg <- sprintf(
  "./figures/jpg/background color %skm %squota %s equal-area latitude bins per-cell all-cells OLS time series.jpg",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

ts_pdf <- sprintf(
  "./figures/pdf/background color %skm %squota %s equal-area latitude bins per-cell all-cells OLS time series.pdf",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

ggsave(ts_jpg, slope_vTime_plot, width = 8, height = 7.4, dpi = 900, device = "jpeg")
ggsave(ts_pdf, slope_vTime_plot, width = 8, height = 7.4, dpi = 300, device = "pdf")

# -----------------------------------------------------------------------
# 8. Boxplot by climate state
# -----------------------------------------------------------------------

slope_data_filtered <- slope_data %>%
  drop_na(slope_value, climate_state) %>%
  filter(climate_state %in% climate_levels) %>%
  mutate(
    climate_state = factor(climate_state, levels = climate_levels),
    slope_type = factor(slope_type, levels = c("Northern", "Southern"))
  )

sample_counts <- slope_data_filtered %>%
  group_by(climate_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    state_with_n = paste0(climate_state, "\n(n=", count, ")"),
    climate_state = factor(climate_state, levels = climate_levels)
  ) %>%
  arrange(climate_state)

slope_data_filtered <- slope_data_filtered %>%
  left_join(
    sample_counts %>% select(climate_state, state_with_n),
    by = "climate_state"
  ) %>%
  mutate(
    state_with_n = factor(state_with_n, levels = sample_counts$state_with_n)
  )

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
  filter(is_outlier)

y_min_box <- min(slope_data_filtered$slope_value, na.rm = TRUE)
y_max_box <- max(slope_data_filtered$slope_value, na.rm = TRUE)

boxplot <- ggplot(
  slope_data_filtered,
  aes(state_with_n, slope_value, fill = slope_type)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -0.1, ymax = 0.1,
    fill = "lightblue", alpha = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(0.75)
  ) +
  # geom_point(
  #   shape = 21,
  #   size = 1,
  #   alpha = 0.6,
  #   color = "black",
  #   position = position_jitterdodge(
  #     jitter.width = 0.2,
  #     dodge.width = 0.75
  #   ),
  #   show.legend = FALSE
  # ) +
  scale_fill_manual(
    values = hemi_cols,
    labels = c("Northern", "Southern"),
    name = "Hemisphere"
  ) +
  labs(
    x = "Climate state",
    y = "Slope value"
  ) +
  coord_cartesian(
    clip = "off",
    xlim = c(1, 5),
    ylim = c(y_min_box, y_max_box)
  ) +
  annotate(
    "text",
    x = 5.81,
    y = (y_max_box + 0.1) / 2,
    label = "Non-modern-type",
    size = 4.5,
    angle = 270
  ) +
  annotate(
    "text",
    x = 5.81,
    y = (y_min_box - 0.1) / 2,
    label = "Modern-type",
    size = 4.5,
    angle = 270
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    ),
    legend.key = element_rect(fill = "white"),
    plot.margin = margin(10, 30, 10, 10)
  )

print(boxplot)

box_jpg <- sprintf(
  "./figures/jpg/%skm %squota %s equal-area latitude bins per-cell all-cells OLS boxplot.jpg",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

box_pdf <- sprintf(
  "./figures/pdf/%skm %squota %s equal-area latitude bins per-cell all-cells OLS boxplot.pdf",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

ggsave(box_jpg, boxplot, width = 7, height = 5, dpi = 300)
ggsave(box_pdf, boxplot, width = 7, height = 5, dpi = 300)

# -----------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------

slope_summary <- slope_cli_df_filter %>%
  group_by(hemisphere) %>%
  summarise(
    good_count = n(),
    greater_than_0 = sum(slope > 0, na.rm = TRUE),
    less_than_0 = sum(slope < 0, na.rm = TRUE),
    equal_to_0 = sum(slope == 0, na.rm = TRUE),
    reverse_LDG = greater_than_0 / good_count * 100,
    normal_LDG = less_than_0 / good_count * 100,
    flat_LDG = equal_to_0 / good_count * 100,
    median_slope = median(slope, na.rm = TRUE),
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE),
    .groups = "drop"
  )

summary_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins per-cell all-cells OLS slope summary.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

write.csv(slope_summary, summary_path, row.names = FALSE)

print(slope_summary)
print(slope_path)
print(slope_cli_path)
print(wil_path)
