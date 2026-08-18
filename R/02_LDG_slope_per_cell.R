# -----------------------------------------------------------------------
# Project: LDG_climate_state
# File: 02_LDG_slope_per_cell.R
# Purpose: Baseline per-cell balanced resampling OLS LDG slopes
#          Baseline = occurrence5_k1_tropical_temperate
# -----------------------------------------------------------------------

# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(palaeoverse)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

# -----------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------

set.seed(123)

n_resamples <- 100  # use 100 for final

baseline_qc <- "occurrence5_k1_tropical_temperate"

baseline_filter <- tibble(
  qc_name = baseline_qc,
  slope_filter_id = "occurrence5",
  qc_type = "occurrence",
  occurrence_min = 5,
  collection_min = NA_real_,
  min_cells_per_latbin = 1,
  zone_requirement = "tropical_temperate",
  require_tropical = TRUE,
  require_temperate = TRUE,
  require_polar = FALSE,
  require_adjacent_tt = TRUE,
  require_min_lat_bins = FALSE
)

percentile_tag <- ifelse(
  length(rich_params$percentiles) > 1,
  "allq",
  rich_params$percentiles[1]
)

analysis_tag <- sprintf(
  "%skm_%squota_%slat_%s",
  params$spacing,
  params$level,
  rich_params$n_lat_bins,
  percentile_tag
)

out_dir <- file.path(
  "./results/LDG_baseline",
  paste0("latbins_", rich_params$n_lat_bins),
  "per_cell_balanced"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 1. Read data
# -----------------------------------------------------------------------

rich_raw <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,
  params$level
)) %>%
  select(-any_of(c("qD_normalized", "qD_raw", "stage_max_qD")))

time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(n = rich_params$n_lat_bins) %>%
  arrange(min)

stopifnot(
  all(c("bin_midpoint", "cell_lat", "collection_no", "t", "nT", "qD") %in% names(rich_raw))
)

lat_zone_lookup <- lat_bins %>%
  mutate(
    abs_lat_bin_mid = round(abs(mid), 6),
    lat_zone = case_when(
      abs_lat_bin_mid < 30 ~ "tropical",
      abs_lat_bin_mid < 60 ~ "temperate",
      abs_lat_bin_mid <= 90 ~ "polar",
      TRUE ~ NA_character_
    )
  ) %>%
  select(bin, abs_lat_bin_mid, lat_zone)

eligible_abs_lat_bins <-
  sort(unique(lat_zone_lookup$abs_lat_bin_mid))

all_stage_hemi <- expand.grid(
  bin_midpoint = sort(unique(rich_raw$bin_midpoint[rich_raw$bin_midpoint <= 486.8500])),
  hemisphere = c("Northern", "Southern"),
  stringsAsFactors = FALSE
) %>%
  mutate(stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)])

# -----------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

fit_cell_ols_safe <- function(df, y_col = "qD_normalized", x_col = "abs_lat") {
  
  df <- df %>%
    filter(!is.na(.data[[y_col]]), !is.na(.data[[x_col]]))
  
  empty <- tibble(
    slope = NA_real_,
    intercept = NA_real_,
    slope_se = NA_real_,
    slope_lower_95 = NA_real_,
    slope_upper_95 = NA_real_,
    p_value = NA_real_,
    r_squared = NA_real_,
    n_cells_model = nrow(df)
  )
  
  if (nrow(df) < 3 || n_distinct(df[[x_col]]) < 2) return(empty)
  
  mod <- tryCatch(
    lm(as.formula(sprintf("%s ~ %s", y_col, x_col)), data = df),
    error = function(e) NULL
  )
  
  if (is.null(mod)) return(empty)
  
  sm <- summary(mod)
  cf <- sm$coefficients
  
  if (!(x_col %in% rownames(cf))) return(empty)
  
  slope <- as.numeric(cf[x_col, "Estimate"])
  se <- as.numeric(cf[x_col, "Std. Error"])
  tibble(
    slope = slope,
    intercept = as.numeric(cf["(Intercept)", "Estimate"]),
    slope_se = se,
    slope_lower_95 = slope - 1.96*se,
    slope_upper_95 = slope + 1.96*se,
    p_value = as.numeric(cf[x_col, "Pr(>|t|)"]),
    r_squared = as.numeric(sm$r.squared),
    n_cells_model = nrow(df)
  )
}

run_balanced_resampling_ols <- function(df) {
  
  df <- df %>%
    filter(!is.na(qD_normalized), !is.na(abs_lat), !is.na(abs_lat_bin_mid))
  
  if (
    nrow(df) < 3 ||
    n_distinct(df$abs_lat) < 2 ||
    n_distinct(df$abs_lat_bin_mid) < 2
  ) {
    return(tibble(
      slope = NA_real_,
      intercept = NA_real_,
      slope_se = NA_real_,
      slope_lower_95 = NA_real_,
      slope_upper_95 = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_cells_model = nrow(df),
      resample_id = seq_len(n_resamples),
      k_sample = NA_real_
    ))
  }
  
  k <- df %>%
    count(abs_lat_bin_mid, name = "n") %>%
    pull(n) %>%
    min(na.rm = TRUE)
  
  split_df <- split(df, df$abs_lat_bin_mid)
  
  map_dfr(seq_len(n_resamples), function(i) {
    
    sampled <- bind_rows(lapply(split_df, function(z) {
      z[
        sample(seq_len(nrow(z)), size = k, replace = FALSE),
        ,
        drop = FALSE
      ]
    }))
    
    fit_cell_ols_safe(sampled) %>%
      mutate(resample_id = i, k_sample = k)
  })
}

summarise_resampled <- function(df) {
  
  s <- df$slope[!is.na(df$slope)]
  
  if (!length(s)) {
    return(tibble(
      slope = NA_real_,
      slope_mean = NA_real_,
      slope_sd = NA_real_,
      slope_lower_95 = NA_real_,
      slope_upper_95 = NA_real_,
      intercept = NA_real_,
      slope_se = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_cells_model = NA_real_,
      n_resamples_valid = 0,
      prop_negative = NA_real_,
      prop_positive = NA_real_,
      k_sample = NA_real_
    ))
  }
  
  tibble(
    slope = median(s),
    slope_mean = mean(s),
    slope_sd = ifelse(length(s) >= 2, sd(s), NA_real_),
    slope_lower_95 = quantile(s, 0.025, names = FALSE),
    slope_upper_95 = quantile(s, 0.975, names = FALSE),
    intercept = median(df$intercept, na.rm = TRUE),
    slope_se = median(df$slope_se, na.rm = TRUE),
    p_value = median(df$p_value, na.rm = TRUE),
    r_squared = median(df$r_squared, na.rm = TRUE),
    n_cells_model = median(df$n_cells_model, na.rm = TRUE),
    n_resamples_valid = length(s),
    prop_negative = mean(s < 0),
    prop_positive = mean(s > 0),
    k_sample = min(df$k_sample, na.rm = TRUE)
  ) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA_real_, .)))
}

filter_valid_latbins_for_qc <- function(rich_df, min_cells_per_latbin) {
  
  if (!nrow(rich_df)) return(rich_df)
  
  valid_latbins <- rich_df %>%
    count(bin_midpoint, hemisphere, abs_lat_bin_mid, name = "n_cells_latbin") %>%
    filter(n_cells_latbin >= min_cells_per_latbin)
  
  if (!nrow(valid_latbins)) return(rich_df[0, ])
  
  rich_df %>%
    inner_join(
      valid_latbins %>% select(bin_midpoint, hemisphere, abs_lat_bin_mid),
      by = c("bin_midpoint", "hemisphere", "abs_lat_bin_mid")
    )
}

empty_hemi_qc <- function() {
  tibble(
    bin_midpoint = numeric(),
    stage = character(),
    hemisphere = character(),
    n_cells = numeric(),
    n_valid_lat_bins = numeric(),
    min_abs_lat_bin_mid = numeric(),
    max_abs_lat_bin_mid = numeric(),
    lat_range_diagnostic = numeric(),
    has_tropical = logical(),
    has_temperate = logical(),
    has_polar = logical(),
    has_adjacent_tt = logical(),
    sum_nT = numeric(),
    mean_nT = numeric(),
    median_nT = numeric(),
    min_nT = numeric(),
    max_nT = numeric(),
    k_min = numeric(),
    k_max = numeric(),
    k_mean = numeric(),
    k_median = numeric(),
    k_cv = numeric(),
    k_evenness = numeric(),
    dominant_lat_bin_prop = numeric(),
    pass_tropical = logical(),
    pass_temperate = logical(),
    pass_polar = logical(),
    pass_adjacent_tt = logical(),
    pass_n_lat_bins = logical(),
    label = character(),
    color = character()
  )
}

classify_stage_hemisphere_qc <- function(rich_for_qc, qc_row, rich_for_sampling = rich_for_qc) {
  
  if (!nrow(rich_for_qc)) return(empty_hemi_qc())
  
  adj <- has_adjacent_bins(rich_for_qc, lat_bins) %>%
    distinct(bin_midpoint, hemisphere, label) %>%
    transmute(bin_midpoint, hemisphere, has_adjacent_tt = label == "good")
  
  hemi <- rich_for_qc %>%
    group_by(bin_midpoint, stage, hemisphere) %>%
    summarise(
      n_cells = n(),
      n_valid_lat_bins = n_distinct(abs_lat_bin_mid),
      min_abs_lat_bin_mid = min(abs_lat_bin_mid, na.rm = TRUE),
      max_abs_lat_bin_mid = max(abs_lat_bin_mid, na.rm = TRUE),
      lat_range_diagnostic = max_abs_lat_bin_mid - min_abs_lat_bin_mid,
      has_tropical = any(lat_zone == "tropical"),
      has_temperate = any(lat_zone == "temperate"),
      has_polar = any(lat_zone == "polar"),
      sum_nT = sum(nT, na.rm = TRUE),
      mean_nT = mean(nT, na.rm = TRUE),
      median_nT = median(nT, na.rm = TRUE),
      min_nT = min(nT, na.rm = TRUE),
      max_nT = max(nT, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(adj, by = c("bin_midpoint", "hemisphere"))
  
  latstats <- rich_for_sampling %>%
    count(bin_midpoint, hemisphere, abs_lat_bin_mid, name = "n_cells_latbin") %>%
    group_by(bin_midpoint, hemisphere) %>%
    complete(
      abs_lat_bin_mid = eligible_abs_lat_bins,
      fill = list(n_cells_latbin = 0L)
    )%>%
    summarise(
      k_min = min(n_cells_latbin),
      k_max = max(n_cells_latbin),
      k_mean = mean(n_cells_latbin),
      k_median = median(n_cells_latbin),
      k_cv = ifelse(k_mean > 0, sd(n_cells_latbin) / k_mean, NA_real_),
      k_evenness = ifelse(k_mean > 0, k_min / k_mean, NA_real_),
      dominant_lat_bin_prop = k_max / sum(n_cells_latbin),
      .groups = "drop"
    )
  
  hemi %>%
    left_join(latstats, by = c("bin_midpoint", "hemisphere")) %>%
    mutate(
      has_tropical = coalesce(has_tropical, FALSE),
      has_temperate = coalesce(has_temperate, FALSE),
      has_polar = coalesce(has_polar, FALSE),
      has_adjacent_tt = coalesce(has_adjacent_tt, FALSE),
      pass_tropical = if (qc_row$require_tropical[[1]]) has_tropical else TRUE,
      pass_temperate = if (qc_row$require_temperate[[1]]) has_temperate else TRUE,
      pass_polar = if (qc_row$require_polar[[1]]) has_polar else TRUE,
      pass_adjacent_tt = if (qc_row$require_adjacent_tt[[1]]) has_adjacent_tt else TRUE,
      pass_n_lat_bins = TRUE,
      label = ifelse(
        pass_tropical & pass_temperate & pass_polar & pass_adjacent_tt,
        "good",
        "bad"
      ),
      color = ifelse(label == "bad", "Bad hemisphere", hemisphere)
    )
}

complete_hemi_qc <- function(hemi_qc, qc_row) {
  
  all_stage_hemi %>%
    left_join(hemi_qc, by = c("bin_midpoint", "stage", "hemisphere")) %>%
    mutate(
      slope_filter_id = qc_row$slope_filter_id[[1]],
      qc_name = qc_row$qc_name[[1]],
      qc_type = qc_row$qc_type[[1]],
      occurrence_min = qc_row$occurrence_min[[1]],
      collection_min = qc_row$collection_min[[1]],
      min_cells_per_latbin = qc_row$min_cells_per_latbin[[1]],
      zone_requirement = qc_row$zone_requirement[[1]],
      require_tropical = qc_row$require_tropical[[1]],
      require_temperate = qc_row$require_temperate[[1]],
      require_polar = qc_row$require_polar[[1]],
      require_adjacent_tt = qc_row$require_adjacent_tt[[1]],
      require_min_lat_bins = qc_row$require_min_lat_bins[[1]],
      label = ifelse(is.na(label), "bad", label),
      color = ifelse(is.na(color), "Bad hemisphere", color)
    )
}

# -----------------------------------------------------------------------
# 3. Prepare baseline richness data
# -----------------------------------------------------------------------

rich_for_slope <- rich_raw %>%
  filter(
    t <= 2 * nT,
    bin_midpoint <= 486.8500,
    nT >= baseline_filter$occurrence_min[[1]]
  ) %>%
  mutate(
    qD_raw = as.numeric(qD),
    stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)],
    bin_index = findInterval(cell_lat, vec = c(lat_bins$min, Inf)),
    bin = lat_bins$bin[bin_index],
    abs_lat = abs(cell_lat),
    hemisphere = ifelse(cell_lat >= 0, "Northern", "Southern")
  ) %>%
  filter(!is.na(bin), !is.na(hemisphere), !is.na(qD_raw)) %>%
  left_join(lat_zone_lookup, by = "bin") %>%
  filter(!is.na(abs_lat_bin_mid), !is.na(lat_zone)) %>%
  group_by(bin_midpoint) %>%
  mutate(qD_normalized = qD_raw * 100 / max(qD_raw, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(qD_normalized)) %>%
  mutate(
    slope_filter_id = baseline_filter$slope_filter_id[[1]],
    qc_name = baseline_filter$qc_name[[1]],
    data_role = "slope"
  )

# -----------------------------------------------------------------------
# 4. Baseline QC labels
# -----------------------------------------------------------------------

rich_for_qc <- filter_valid_latbins_for_qc(
  rich_for_slope,
  baseline_filter$min_cells_per_latbin[[1]]
) %>%
  mutate(
    qc_name = baseline_filter$qc_name[[1]],
    slope_filter_id = baseline_filter$slope_filter_id[[1]],
    data_role = "QC"
  )

hemi_QC_summary <- classify_stage_hemisphere_qc(
  rich_for_qc,
  baseline_filter,
  rich_for_slope
) %>%
  complete_hemi_qc(baseline_filter)

# -----------------------------------------------------------------------
# 5. Sampling-profile dissimilarity
# -----------------------------------------------------------------------

sampling_profile_dissim <- rich_for_slope %>%
  count(
    slope_filter_id,
    bin_midpoint,
    hemisphere,
    abs_lat_bin_mid,
    name = "n_cells_latbin"
  ) %>%
  group_by(slope_filter_id, bin_midpoint) %>%
  complete(
    hemisphere = c("Northern", "Southern"),
    abs_lat_bin_mid,
    fill = list(n_cells_latbin = 0)
  ) %>%
  ungroup() %>%
  group_by(slope_filter_id, bin_midpoint, hemisphere) %>%
  mutate(
    p_cells = ifelse(
      sum(n_cells_latbin, na.rm = TRUE) > 0,
      n_cells_latbin / sum(n_cells_latbin, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  select(
    slope_filter_id,
    bin_midpoint,
    hemisphere,
    abs_lat_bin_mid,
    p_cells
  ) %>%
  pivot_wider(
    names_from = hemisphere,
    values_from = p_cells,
    names_prefix = "p_"
  ) %>%
  group_by(slope_filter_id, bin_midpoint) %>%
  summarise(
    sampling_profile_dissim = 0.5 * sum(
      abs(coalesce(p_Northern, 0) - coalesce(p_Southern, 0)),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# -----------------------------------------------------------------------
# 6. Per-cell balanced resampling OLS slopes
# -----------------------------------------------------------------------
set.seed(123)
percell_raw_resamples <- rich_for_slope %>%
  group_by(bin_midpoint, stage, hemisphere) %>%
  group_modify(~ run_balanced_resampling_ols(.x)) %>%
  ungroup() %>%
  mutate(slope_filter_id = baseline_filter$slope_filter_id[[1]])

LDG_slope <- percell_raw_resamples %>%
  group_by(bin_midpoint, stage, hemisphere) %>%
  group_modify(~ summarise_resampled(.x)) %>%
  ungroup() %>%
  right_join(all_stage_hemi, by = c("bin_midpoint", "stage", "hemisphere")) %>%
  mutate(
    slope_filter_id = baseline_filter$slope_filter_id[[1]],
    slope_qc_type = baseline_filter$qc_type[[1]],
    slope_occurrence_min = baseline_filter$occurrence_min[[1]],
    slope_collection_min = baseline_filter$collection_min[[1]],
    method_group = "per_cell_balanced_resampling_OLS",
    slope_metric = "median_resampled_slope",
    quantile = NA_character_
  ) %>%
  left_join(
    hemi_QC_summary,
    by = c("slope_filter_id", "bin_midpoint", "stage", "hemisphere")
  ) %>%
  left_join(
    sampling_profile_dissim,
    by = c("slope_filter_id", "bin_midpoint")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label),
    color = ifelse(is.na(color), "Bad hemisphere", color),
    slope_90deg = slope * 90,
    slope_sampled_range = slope * lat_range_diagnostic,
    interval_width = slope_upper_95 - slope_lower_95,
    slope_direction = case_when(
      slope < 0 ~ "normal_LDG",
      slope > 0 ~ "reverse_LDG",
      slope == 0 ~ "flat_LDG",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(hemisphere, bin_midpoint)

# -----------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------

# Keep the same output names as before for downstream plotting scripts
out_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

raw_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS raw resamples.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

# Add old-compatible aliases before saving
LDG_slope <- LDG_slope %>%
  mutate(
    method = "grid_cell_balanced_resampling_OLS",
    lat_range = lat_range_diagnostic,
    min_abs_lat = min_abs_lat_bin_mid,
    max_abs_lat = max_abs_lat_bin_mid
  )

write.csv(
  LDG_slope,
  out_path,
  row.names = FALSE
)

write.csv(
  percell_raw_resamples,
  raw_path,
  row.names = FALSE
)

# -----------------------------------------------------------------------
# 8. Quick check
# -----------------------------------------------------------------------

print(
  LDG_slope %>%
    group_by(hemisphere, label) %>%
    summarise(
      n = n(),
      n_valid = sum(!is.na(slope)),
      median_slope = median(slope, na.rm = TRUE),
      median_interval_width = median(interval_width, na.rm = TRUE),
      median_k_sample = median(k_sample, na.rm = TRUE),
      median_k_median = median(k_median, na.rm = TRUE),
      median_k_cv = median(k_cv, na.rm = TRUE),
      .groups = "drop"
    )
)

print(out_path)
print(raw_path)

message("Baseline per-cell balanced resampling OLS finished.")