# -----------------------------------------------------------------------
# Project: LDG_climate_state
# File: 05_high_latitude_coverage_summary_FigS3.R
# Purpose: Summarise high-latitude coverage for LDG slope estimates
# -----------------------------------------------------------------------

# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(palaeoverse)
  library(deeptime)
  library(patchwork)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

# -----------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------

occurrence_min <- 5

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
dir.create("./figures/jpg", recursive = TRUE, showWarnings = FALSE)
dir.create("./figures/pdf", recursive = TRUE, showWarnings = FALSE)

hemi_cols <- c("Northern" = "#0072B2", "Southern" = "#E69F00")

# -----------------------------------------------------------------------
# 2. Read and prepare richness data
# -----------------------------------------------------------------------

rich_df <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing, params$level
))

time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(n = rich_params$n_lat_bins) %>%
  arrange(min)

lat_zone_lookup <- lat_bins %>%
  mutate(
    lat_bin_mid = mid,
    abs_lat_bin_mid = round(abs(mid), 6),
    lat_zone = case_when(
      abs_lat_bin_mid < 30 ~ "tropical",
      abs_lat_bin_mid < 60 ~ "temperate",
      abs_lat_bin_mid <= 90 ~ "polar",
      TRUE ~ NA_character_
    )
  ) %>%
  select(bin, lat_bin_mid, abs_lat_bin_mid, lat_zone)

rich_df <- rich_df %>%
  filter(bin_midpoint <= 486.8500) %>%
  mutate(
    stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)],
    completeness = ifelse(nT >= occurrence_min & t <= 2 * nT, "Complete", "Incomplete"),
    bin_index = findInterval(cell_lat, vec = c(lat_bins$min, Inf)),
    bin = lat_bins$bin[bin_index],
    abs_lat = abs(cell_lat),
    hemisphere = case_when(
      cell_lat >= 0 ~ "Northern",
      cell_lat < 0 ~ "Southern",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(lat_zone_lookup, by = "bin") %>%
  filter(
    !is.na(bin),
    !is.na(abs_lat),
    !is.na(abs_lat_bin_mid),
    !is.na(lat_zone),
    !is.na(hemisphere)
  )

# -----------------------------------------------------------------------
# 3. Baseline QC: occurrence5_k1_tropical_temperate
# -----------------------------------------------------------------------

rich_complete <- rich_df %>%
  filter(completeness == "Complete")

adjacent_df <- has_adjacent_bins(rich_complete, lat_bins) %>%
  distinct(bin_midpoint, hemisphere, label) %>%
  transmute(
    bin_midpoint,
    hemisphere,
    has_adjacent_tt = label == "good"
  )

qc_label_df <- rich_complete %>%
  group_by(bin_midpoint, hemisphere) %>%
  summarise(
    has_tropical = any(lat_zone == "tropical", na.rm = TRUE),
    has_temperate = any(lat_zone == "temperate", na.rm = TRUE),
    has_polar = any(lat_zone == "polar", na.rm = TRUE),
    n_valid_lat_bins = n_distinct(abs_lat_bin_mid),
    .groups = "drop"
  ) %>%
  left_join(adjacent_df, by = c("bin_midpoint", "hemisphere")) %>%
  mutate(
    has_adjacent_tt = coalesce(has_adjacent_tt, FALSE),
    label = ifelse(
      has_tropical & has_temperate & has_adjacent_tt,
      "good",
      "bad"
    )
  )

rich_df <- rich_df %>%
  left_join(
    qc_label_df %>%
      select(
        bin_midpoint,
        hemisphere,
        label,
        has_tropical,
        has_temperate,
        has_polar,
        has_adjacent_tt,
        n_valid_lat_bins
      ),
    by = c("bin_midpoint", "hemisphere")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label),
    label = factor(label, levels = c("good", "bad"))
  )

# -----------------------------------------------------------------------
# 4. Stage × hemisphere high-latitude coverage
# -----------------------------------------------------------------------

coverage_df <- rich_df %>%
  group_by(bin_midpoint, stage, hemisphere, label) %>%
  summarise(
    n_cells = n(),
    n_lat_bins = n_distinct(abs_lat_bin_mid),
    min_abs_lat = min(abs_lat, na.rm = TRUE),
    max_abs_lat = max(abs_lat, na.rm = TRUE),
    min_abs_lat_bin_mid = min(abs_lat_bin_mid, na.rm = TRUE),
    max_abs_lat_bin_mid = max(abs_lat_bin_mid, na.rm = TRUE),
    n_cells_ge_60 = sum(abs_lat >= 60, na.rm = TRUE),
    n_cells_ge_66_5 = sum(abs_lat >= 66.5, na.rm = TRUE),
    n_cells_ge_70 = sum(abs_lat >= 70, na.rm = TRUE),
    n_cells_ge_75 = sum(abs_lat >= 75, na.rm = TRUE),
    prop_cells_ge_60 = n_cells_ge_60 / n_cells,
    prop_cells_ge_66_5 = n_cells_ge_66_5 / n_cells,
    prop_cells_ge_70 = n_cells_ge_70 / n_cells,
    prop_cells_ge_75 = n_cells_ge_75 / n_cells,
    reaches_60 = max_abs_lat >= 60,
    reaches_66_5 = max_abs_lat >= 66.5,
    reaches_70 = max_abs_lat >= 70,
    reaches_75 = max_abs_lat >= 75,
    .groups = "drop"
  )

coverage_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins high-latitude coverage by stage-hemisphere.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(coverage_df, coverage_path, row.names = FALSE)

# -----------------------------------------------------------------------
# 5. Summary: all good stage × hemisphere datasets
# -----------------------------------------------------------------------

coverage_summary_good <- coverage_df %>%
  filter(label == "good") %>%
  summarise(
    n_good_stage_hemisphere = n(),
    median_max_abs_lat = median(max_abs_lat, na.rm = TRUE),
    mean_max_abs_lat = mean(max_abs_lat, na.rm = TRUE),
    q25_max_abs_lat = quantile(max_abs_lat, 0.25, na.rm = TRUE),
    q75_max_abs_lat = quantile(max_abs_lat, 0.75, na.rm = TRUE),
    n_reaches_60 = sum(reaches_60, na.rm = TRUE),
    n_reaches_66_5 = sum(reaches_66_5, na.rm = TRUE),
    n_reaches_70 = sum(reaches_70, na.rm = TRUE),
    n_reaches_75 = sum(reaches_75, na.rm = TRUE),
    prop_reaches_60 = n_reaches_60 / n_good_stage_hemisphere,
    prop_reaches_66_5 = n_reaches_66_5 / n_good_stage_hemisphere,
    prop_reaches_70 = n_reaches_70 / n_good_stage_hemisphere,
    prop_reaches_75 = n_reaches_75 / n_good_stage_hemisphere,
    median_prop_cells_ge_60 = median(prop_cells_ge_60, na.rm = TRUE),
    median_prop_cells_ge_66_5 = median(prop_cells_ge_66_5, na.rm = TRUE),
    median_prop_cells_ge_70 = median(prop_cells_ge_70, na.rm = TRUE),
    median_prop_cells_ge_75 = median(prop_cells_ge_75, na.rm = TRUE)
  )

summary_good_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins high-latitude coverage summary good.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(coverage_summary_good, summary_good_path, row.names = FALSE)
print(coverage_summary_good)

# -----------------------------------------------------------------------
# 6. Summary by hemisphere
# -----------------------------------------------------------------------

coverage_summary_hemi <- coverage_df %>%
  filter(label == "good") %>%
  group_by(hemisphere) %>%
  summarise(
    n_good_stage_hemisphere = n(),
    n_reaches_60 = sum(reaches_60, na.rm = TRUE),
    n_reaches_66_5 = sum(reaches_66_5, na.rm = TRUE),
    n_reaches_70 = sum(reaches_70, na.rm = TRUE),
    n_reaches_75 = sum(reaches_75, na.rm = TRUE),
    prop_reaches_60 = n_reaches_60 / n_good_stage_hemisphere,
    prop_reaches_66_5 = n_reaches_66_5 / n_good_stage_hemisphere,
    prop_reaches_70 = n_reaches_70 / n_good_stage_hemisphere,
    prop_reaches_75 = n_reaches_75 / n_good_stage_hemisphere,
    median_max_abs_lat = median(max_abs_lat, na.rm = TRUE),
    mean_max_abs_lat = mean(max_abs_lat, na.rm = TRUE),
    q25_max_abs_lat = quantile(max_abs_lat, 0.25, na.rm = TRUE),
    q75_max_abs_lat = quantile(max_abs_lat, 0.75, na.rm = TRUE),
    median_prop_cells_ge_60 = median(prop_cells_ge_60, na.rm = TRUE),
    median_prop_cells_ge_66_5 = median(prop_cells_ge_66_5, na.rm = TRUE),
    median_prop_cells_ge_70 = median(prop_cells_ge_70, na.rm = TRUE),
    median_prop_cells_ge_75 = median(prop_cells_ge_75, na.rm = TRUE),
    .groups = "drop"
  )

summary_hemi_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins high-latitude coverage summary by hemisphere.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(coverage_summary_hemi, summary_hemi_path, row.names = FALSE)
print(coverage_summary_hemi)

# -----------------------------------------------------------------------
# 7. Paired NH-SH high-latitude coverage summary
# -----------------------------------------------------------------------

coverage_pair <- coverage_df %>%
  filter(label == "good") %>%
  select(
    bin_midpoint, stage, hemisphere,
    max_abs_lat,
    prop_cells_ge_60,
    prop_cells_ge_66_5,
    prop_cells_ge_70,
    prop_cells_ge_75,
    reaches_60,
    reaches_66_5,
    reaches_70,
    reaches_75
  ) %>%
  pivot_wider(
    names_from = hemisphere,
    values_from = c(
      max_abs_lat,
      prop_cells_ge_60,
      prop_cells_ge_66_5,
      prop_cells_ge_70,
      prop_cells_ge_75,
      reaches_60,
      reaches_66_5,
      reaches_70,
      reaches_75
    ),
    names_sep = "_"
  ) %>%
  filter(!is.na(max_abs_lat_Northern), !is.na(max_abs_lat_Southern)) %>%
  mutate(
    both_reach_60 = reaches_60_Northern & reaches_60_Southern,
    both_reach_66_5 = reaches_66_5_Northern & reaches_66_5_Southern,
    both_reach_70 = reaches_70_Northern & reaches_70_Southern,
    both_reach_75 = reaches_75_Northern & reaches_75_Southern,
    max_abs_lat_diff = abs(max_abs_lat_Northern - max_abs_lat_Southern),
    prop_cells_ge_60_diff = abs(prop_cells_ge_60_Northern - prop_cells_ge_60_Southern),
    prop_cells_ge_70_diff = abs(prop_cells_ge_70_Northern - prop_cells_ge_70_Southern)
  )

coverage_pair_summary <- coverage_pair %>%
  summarise(
    n_paired_good_stages = n(),
    prop_both_reach_60 = mean(both_reach_60, na.rm = TRUE),
    prop_both_reach_66_5 = mean(both_reach_66_5, na.rm = TRUE),
    prop_both_reach_70 = mean(both_reach_70, na.rm = TRUE),
    prop_both_reach_75 = mean(both_reach_75, na.rm = TRUE),
    median_max_abs_lat_diff = median(max_abs_lat_diff, na.rm = TRUE),
    mean_max_abs_lat_diff = mean(max_abs_lat_diff, na.rm = TRUE),
    median_prop_cells_ge_60_diff = median(prop_cells_ge_60_diff, na.rm = TRUE),
    median_prop_cells_ge_70_diff = median(prop_cells_ge_70_diff, na.rm = TRUE)
  )

pair_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins paired high-latitude coverage.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

pair_summary_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins paired high-latitude coverage summary.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(coverage_pair, pair_path, row.names = FALSE)
write.csv(coverage_pair_summary, pair_summary_path, row.names = FALSE)
print(coverage_pair_summary)

# -----------------------------------------------------------------------
# 8. Plot maximum palaeolatitudinal coverage through time
# -----------------------------------------------------------------------

data(periods)
data(epochs)

major_boundaries <- periods$max_age
x_max_val <- max(time_bins$max_ma, na.rm = TRUE)
y_breaks <- seq(0, 90, 15)

cov_plot <- coverage_df %>%
  mutate(
    label = factor(label, levels = c("good", "bad")),
    point_fill = ifelse(label == "good", hemisphere, "Bad")
  )

make_sum_txt <- function(hemi_use) {
  coverage_summary_hemi %>%
    filter(hemisphere == hemi_use) %>%
    transmute(
      txt = paste0(
        "≥60°: ", n_reaches_60, "/", n_good_stage_hemisphere,
        " (", round(prop_reaches_60 * 100, 1), "%)\n",
        "≥70°: ", n_reaches_70, "/", n_good_stage_hemisphere,
        " (", round(prop_reaches_70 * 100, 1), "%)\n",
        "Median max: ", round(median_max_abs_lat, 1), "°"
      )
    ) %>%
    pull(txt)
}

draw_cov_panel <- function(hemi_use, hemi_col, tag_lab, show_x = FALSE) {
  
  df <- cov_plot %>% filter(hemisphere == hemi_use)
  
  p <- ggplot(df, aes(x = bin_midpoint, y = max_abs_lat)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = 60, ymax = 90,
      fill = "grey80", alpha = 0.35
    ) +
    geom_vline(
      xintercept = major_boundaries,
      color = "black",
      linewidth = 0.4,
      alpha = 0.8
    ) +
    geom_hline(yintercept = 60, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
    geom_hline(yintercept = 70, linetype = "dotted", colour = "grey25", linewidth = 0.45) +
    geom_line(linewidth = 0.9, colour = hemi_col, alpha = 0.8, na.rm = TRUE) +
    geom_point(
      data = df %>% filter(label == "good"),
      shape = 21,
      size = 2.3,
      stroke = 0.6,
      colour = "black",
      fill = hemi_col
    ) +
    geom_point(
      data = df %>% filter(label == "bad"),
      shape = 21,
      size = 2.3,
      stroke = 0.6,
      colour = hemi_col,
      fill = "white"
    ) +
    annotate(
      "label",
      x = 20,
      y = 35,
      label = make_sum_txt(hemi_use),
      hjust = 1,
      vjust = 1,
      size = 3.0,
      linewidth = 0.3,
      fill = "white",
      colour = "black"
    ) +
    scale_x_reverse(
      name = ifelse(show_x, "Time (Ma)", ""),
      limits = c(x_max_val, 0),
      breaks = seq(500, 0, -50),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 90),
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    labs(
      x = if (show_x) "Time (Ma)" else "",
      y = "Max abs. palaeolat. (°)",
      tag = tag_lab,
      subtitle = paste(hemi_use, "Hemisphere")
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, colour = "black"),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      plot.tag = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 13, face = "bold"),
      plot.margin = margin(5, 8, 2, 8)
    )
  
  if (show_x) {
    p <- p +
      coord_geo(xlim =c(x_max_val,0),
                pos = as.list(rep("bottom", 2)),
                dat = list("periods", "era"),
                height = list(unit(1.35, "lines"), unit(1.35, "lines")),
                #fill="darkgrey",
                lab_color="black",
                rot = list(0, 0), 
                # size = list(4, 4), 
                abbrv =list(TRUE,FALSE)) +
      theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5))
  } else {
    p <- p +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  p
}

P_north_cov <- draw_cov_panel("Northern", hemi_cols["Northern"], "(a)", FALSE)
P_south_cov <- draw_cov_panel("Southern", hemi_cols["Southern"], "(b)", TRUE)

p_cov_final <- P_north_cov / P_south_cov +
  plot_layout(heights = c(1, 1))

print(p_cov_final)

cov_jpg <- sprintf(
  "./figures/jpg/FigS3_%skm %squota %s equal-area latitude bins maximum high-latitude coverage split.jpg",
  params$spacing, params$level, rich_params$n_lat_bins
)

cov_pdf <- sprintf(
  "./figures/pdf/%skm %squota %s equal-area latitude bins maximum high-latitude coverage split.pdf",
  params$spacing, params$level, rich_params$n_lat_bins
)

ggsave(cov_jpg, p_cov_final, width = 8, height = 7, dpi = 300)
ggsave(cov_pdf, p_cov_final, width = 8, height = 7, device = "pdf")

# -----------------------------------------------------------------------
# 9. One-line summary for manuscript / response
# -----------------------------------------------------------------------

summary_text <- coverage_summary_good %>%
  transmute(
    txt = paste0(
      "Among quality-controlled stage–hemisphere datasets (n = ",
      n_good_stage_hemisphere,
      "), ",
      round(prop_reaches_60 * 100, 1), "% reached >=60°, ",
      round(prop_reaches_66_5 * 100, 1), "% reached >=66.5°, ",
      round(prop_reaches_70 * 100, 1), "% reached >=70°, and ",
      round(prop_reaches_75 * 100, 1), "% reached >=75°. ",
      "The median maximum absolute palaeolatitude was ",
      round(median_max_abs_lat, 1), "°."
    )
  ) %>%
  pull(txt)

cat("\n", summary_text, "\n")

print(coverage_path)
print(summary_good_path)
print(summary_hemi_path)
print(pair_path)
print(pair_summary_path)
print(cov_jpg)
print(cov_pdf)