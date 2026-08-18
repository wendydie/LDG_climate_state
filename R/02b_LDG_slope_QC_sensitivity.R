# -----------------------------------------------------------------------
# Project: LDG_climate_state
# File: 02b_LDG_slope_QC_sensitivity.R
# Purpose: Optimised QC sensitivity for percentile, balanced per-cell,
#          and all-cell OLS LDG slopes
# -----------------------------------------------------------------------

# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(palaeoverse)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(deeptime)
  library(ggh4x)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

# -----------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------

set.seed(123)
filter_seeds <- c(
  occurrence5  = 123L,
  occurrence10 = 124L,
  collection5  = 125L,
  collection10 = 126L
)
n_resamples <- 100          # use 100 for final
save_raw_resamples <- FALSE # TRUE only for debugging
min_time_series_retention <- 0.5

baseline_qc <- "occurrence5_k1_tropical_temperate"

climate_levels <- c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse")

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

era_cols <- c(
  "Paleozoic" = "#9BBA7F",
  "Mesozoic"  = "#67C5CA",
  "Cenozoic"  = "#F2D2A2"
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
  "./results/LDG_QC_sensitivity",
  paste0("latbins_", rich_params$n_lat_bins),
  "all_methods_QC_only"
)

fig_dir <- file.path(
  "./figures/LDG_QC_sensitivity",
  paste0("latbins_", rich_params$n_lat_bins),
  "all_methods_QC_only"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "jpg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "time_series_jpg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "time_series_pdf"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 1. Data
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

climate_states <- read.csv("./data/climate_states.csv") %>%
  filter(bottom <= 486.8500, top >= 0) %>%
  mutate(
    climate_state = factor(climate_state, levels = climate_levels),
    climate_color = climate_colors[as.character(climate_state)]
  )

time_bins_use <- time_bins %>%
  filter(min_ma <= 486.8500, max_ma >= 0)

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

eligible_abs_lat_bins <- sort(
  unique(lat_zone_lookup$abs_lat_bin_mid)
)

all_stage_hemi <- expand.grid(
  bin_midpoint = sort(unique(rich_raw$bin_midpoint[rich_raw$bin_midpoint <= 486.8500])),
  hemisphere = c("Northern", "Southern"),
  stringsAsFactors = FALSE
) %>%
  mutate(stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)])

# -----------------------------------------------------------------------
# 2. QC scenarios, slope filters, and method specs
# -----------------------------------------------------------------------
qc_scenarios <- bind_rows(
  crossing(
    qc_type = "occurrence",
    occurrence_min = c(5, 10),
    collection_min = NA_real_,
    min_cells_per_latbin = c(1, 2, 3),
    zone_requirement = c("tropical_temperate", "tropical_temperate_polar")
  ),
  crossing(
    qc_type = "collection",
    occurrence_min = NA_real_,
    collection_min = c(5, 10),
    min_cells_per_latbin = c(1, 2, 3),
    zone_requirement = c("tropical_temperate", "tropical_temperate_polar")
  )
) %>%
  mutate(
    require_tropical = TRUE,
    require_temperate = TRUE,
    require_polar = zone_requirement == "tropical_temperate_polar",
    require_adjacent_tt = TRUE,
    require_min_lat_bins = FALSE,
    slope_filter_id = case_when(
      qc_type == "occurrence" ~ paste0("occurrence", occurrence_min),
      qc_type == "collection" ~ paste0("collection", collection_min)
    ),
    qc_name = case_when(
      qc_type == "occurrence" ~ paste0(
        "occurrence", occurrence_min,
        "_k", min_cells_per_latbin,
        "_", zone_requirement
      ),
      qc_type == "collection" ~ paste0(
        "collection", collection_min,
        "_k", min_cells_per_latbin,
        "_", zone_requirement
      )
    )
  )
slope_filters <- qc_scenarios %>%
  distinct(slope_filter_id, qc_type, occurrence_min, collection_min) %>%
  arrange(qc_type, occurrence_min, collection_min)

method_specs <- bind_rows(
  tibble(method_group = "percentile_latbin_OLS", slope_metric = as.character(rich_params$percentiles)),
  tibble(method_group = "per_cell_balanced_resampling_OLS", slope_metric = "median_resampled_slope"),
  tibble(method_group = "per_cell_all_cells_OLS", slope_metric = "all_cells_slope")
)

method_specs_main <- bind_rows(
  tibble(method_group = "percentile_latbin_OLS", slope_metric = "q75"),
  tibble(method_group = "per_cell_balanced_resampling_OLS", slope_metric = "median_resampled_slope"),
  tibble(method_group = "per_cell_all_cells_OLS", slope_metric = "all_cells_slope")
)

stopifnot("slope_filter_id" %in% names(qc_scenarios))
stopifnot("slope_filter_id" %in% names(slope_filters))

write.csv(qc_scenarios, file.path(out_dir, paste0(analysis_tag, "_QC_scenarios.csv")), row.names = FALSE)
write.csv(slope_filters, file.path(out_dir, paste0(analysis_tag, "_slope_filters.csv")), row.names = FALSE)

# -----------------------------------------------------------------------
# 3. Helpers
# -----------------------------------------------------------------------

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_median <- function(x) {
  x <- safe_numeric(x)
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else median(x)
}

safe_iqr <- function(x) {
  x <- safe_numeric(x)
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else IQR(x)
}

safe_cor <- function(x, y) {
  x <- safe_numeric(x)
  y <- safe_numeric(y)
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

safe_wilcox <- function(x, y, min_group_n = 2) {
  x <- safe_numeric(x)
  y <- safe_numeric(y)
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  if (length(x) < min_group_n || length(y) < min_group_n) {
    return(list(
      p_value = NA_real_,
      w_statistic = NA_real_,
      test_note = paste0("Skipped: n < ", min_group_n)
    ))
  }
  
  wt <- suppressWarnings(wilcox.test(x, y, alternative = "less", exact = FALSE))
  
  list(
    p_value = wt$p.value,
    w_statistic = as.numeric(wt$statistic),
    test_note = "Tested"
  )
}

cor_fun <- function(x, y) {
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 3) {
    return(tibble(n = sum(ok), estimate = NA_real_, p_value = NA_real_))
  }
  
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  
  tibble(n = sum(ok), estimate = unname(z$estimate), p_value = z$p.value)
}

cor_one <- function(data, x, y, test_name) {
  cor_fun(data[[x]], data[[y]]) %>%
    mutate(test = test_name, x_var = x, y_var = y) %>%
    select(test, x_var, y_var, n, estimate, p_value)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

clean_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

safe_range_expand <- function(x, expand_frac = 0.08) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(-1, 1))
  xr <- range(x)
  if (xr[1] == xr[2]) xr + c(-0.5, 0.5) else xr + c(-1, 1) * diff(xr) * expand_frac
}

# -----------------------------------------------------------------------
# 4. Prepare data and QC
# -----------------------------------------------------------------------

prepare_richness_data <- function(rich_df, filter_row) {
  
  out <- rich_df %>%
    filter(t <= 2 * nT, bin_midpoint <= 486.8500)
  
  if (filter_row$qc_type[[1]] == "occurrence") {
    out <- out %>% filter(nT >= filter_row$occurrence_min[[1]])
  }
  
  if (filter_row$qc_type[[1]] == "collection") {
    out <- out %>% filter(collection_no >= filter_row$collection_min[[1]])
  }
  
  if (!nrow(out)) return(out)
  
  out %>%
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
    filter(!is.na(qD_normalized))
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

classify_stage_hemisphere_qc <- function(rich_for_qc, qc_row,rich_for_sampling = rich_for_qc) {
  
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
    ) %>%
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
# 5. Slope functions
# -----------------------------------------------------------------------

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
  
  if (nrow(df) < 3 ||
      n_distinct(df$abs_lat) < 2 ||
      n_distinct(df$abs_lat_bin_mid) < 2) {
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
      z[sample(seq_len(nrow(z)), size = k, replace =  FALSE), , drop = FALSE]
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

calc_percentile_slopes_base <- function(rich_for_slope, slope_filter_row) {
  
  rich_summary <- rich_for_slope %>%
    group_by(bin_midpoint, stage, hemisphere, abs_lat_bin_mid) %>%
    summarise(
      q50 = quantile(qD_normalized, 0.50, na.rm = TRUE),
      q60 = quantile(qD_normalized, 0.60, na.rm = TRUE),
      q75 = quantile(qD_normalized, 0.75, na.rm = TRUE),
      q90 = quantile(qD_normalized, 0.90, na.rm = TRUE),
      q95 = quantile(qD_normalized, 0.95, na.rm = TRUE),
      n_cells_latbin = n(),
      .groups = "drop"
    )
  
  res <- map_dfr(rich_params$percentiles, function(q) {
    rich_summary %>%
      group_by(bin_midpoint, stage, hemisphere) %>%
      group_modify(~ {
        df <- .x %>%
          filter(
            !is.na(.data[[q]]),
            !is.na(abs_lat_bin_mid)
          )
        
        # 至少需要3个纬度带，OLS才能估计残差和置信区间
        if (
          nrow(df) < 3 ||
          n_distinct(df$abs_lat_bin_mid) < 3
        ) {
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
          lm(
            as.formula(sprintf("%s ~ abs_lat_bin_mid", q)),
            data = df
          ),
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
        
        slope_ci <- tryCatch(
          confint(
            mod,
            parm = "abs_lat_bin_mid",
            level = 0.95
          ),
          error = function(e) c(NA_real_, NA_real_)
        )
        
        tibble(
          slope = unname(coef(mod)["abs_lat_bin_mid"]),
          intercept = unname(coef(mod)["(Intercept)"]),
          slope_se = unname(
            cf["abs_lat_bin_mid", "Std. Error"]
          ),
          slope_lower_95 = unname(slope_ci[1]),
          slope_upper_95 = unname(slope_ci[2]),
          p_value = unname(
            cf["abs_lat_bin_mid", "Pr(>|t|)"]
          ),
          r_squared = unname(sm$r.squared),
          n_cells_model = nrow(df)
        )
      }) %>%
      ungroup() %>%
      mutate(
        method_group = "percentile_latbin_OLS",
        slope_metric = q,
        quantile = q
      )
  })
  
  all_grid <- expand.grid(
    bin_midpoint = sort(unique(rich_raw$bin_midpoint[rich_raw$bin_midpoint <= 486.8500])),
    hemisphere = c("Northern", "Southern"),
    slope_metric = as.character(rich_params$percentiles),
    stringsAsFactors = FALSE
  ) %>%
    mutate(stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)])
  
  res %>%
    right_join(
      all_grid,
      by = c("bin_midpoint", "stage", "hemisphere", "slope_metric")
    ) %>%
    mutate(
      slope_filter_id = slope_filter_row$slope_filter_id[[1]],
      slope_qc_type = slope_filter_row$qc_type[[1]],
      slope_occurrence_min = slope_filter_row$occurrence_min[[1]],
      slope_collection_min = slope_filter_row$collection_min[[1]],
      method_group = "percentile_latbin_OLS",
      quantile = slope_metric,
      slope = safe_numeric(slope),
      intercept = safe_numeric(intercept),
      slope_mean = NA_real_,
      slope_sd = NA_real_,
      
      # 保留前面OLS计算得到的结果
      slope_lower_95 = safe_numeric(slope_lower_95),
      slope_upper_95 = safe_numeric(slope_upper_95),
      slope_se = safe_numeric(slope_se),
      p_value = safe_numeric(p_value),
      r_squared = safe_numeric(r_squared),
      n_cells_model = safe_numeric(n_cells_model),
      
      n_resamples_valid = NA_real_,
      prop_negative = NA_real_,
      prop_positive = NA_real_,
      k_sample = NA_real_
    )
}

calc_percell_all_slopes_base <- function(rich_for_slope, slope_filter_row) {
  
  rich_for_slope %>%
    group_by(bin_midpoint, stage, hemisphere) %>%
    group_modify(~ fit_cell_ols_safe(.x)) %>%
    ungroup() %>%
    right_join(all_stage_hemi, by = c("bin_midpoint", "stage", "hemisphere")) %>%
    mutate(
      slope_filter_id = slope_filter_row$slope_filter_id[[1]],
      slope_qc_type = slope_filter_row$qc_type[[1]],
      slope_occurrence_min = slope_filter_row$occurrence_min[[1]],
      slope_collection_min = slope_filter_row$collection_min[[1]],
      method_group = "per_cell_all_cells_OLS",
      slope_metric = "all_cells_slope",
      quantile = NA_character_,
      slope_mean = NA_real_,
      slope_sd = NA_real_,
      n_resamples_valid = NA_real_,
      prop_negative = ifelse(!is.na(slope), as.numeric(slope < 0), NA_real_),
      prop_positive = ifelse(!is.na(slope), as.numeric(slope > 0), NA_real_),
      k_sample = NA_real_
    )
}

calc_percell_balanced_slopes_base <- function(rich_for_slope, slope_filter_row) {
  
  raw <- rich_for_slope %>%
    group_by(bin_midpoint, stage, hemisphere) %>%
    group_modify(~ run_balanced_resampling_ols(.x)) %>%
    ungroup() %>%
    mutate(slope_filter_id = slope_filter_row$slope_filter_id[[1]])
  
  slope <- raw %>%
    group_by(bin_midpoint, stage, hemisphere) %>%
    group_modify(~ summarise_resampled(.x)) %>%
    ungroup() %>%
    right_join(all_stage_hemi, by = c("bin_midpoint", "stage", "hemisphere")) %>%
    mutate(
      slope_filter_id = slope_filter_row$slope_filter_id[[1]],
      slope_qc_type = slope_filter_row$qc_type[[1]],
      slope_occurrence_min = slope_filter_row$occurrence_min[[1]],
      slope_collection_min = slope_filter_row$collection_min[[1]],
      method_group = "per_cell_balanced_resampling_OLS",
      slope_metric = "median_resampled_slope",
      quantile = NA_character_
    )
  
  list(slope = slope, raw = raw)
}

# -----------------------------------------------------------------------
# 6. Run slope filters and QC labels separately
# -----------------------------------------------------------------------

run_one_slope_filter <- function(slope_filter_row) {
  
  slope_filter_row <- as_tibble(slope_filter_row)
  
  filter_id <- slope_filter_row$slope_filter_id[[1]]
  
  rich_for_slope <- prepare_richness_data(
    rich_raw,
    slope_filter_row
  ) %>%
    mutate(
      slope_filter_id = filter_id,
      data_role = "slope"
    )
  
  perc <- calc_percentile_slopes_base(
    rich_for_slope,
    slope_filter_row
  )
  
  allc <- calc_percell_all_slopes_base(
    rich_for_slope,
    slope_filter_row
  )
  
  set.seed(filter_seeds[[filter_id]])
  
  bal <- calc_percell_balanced_slopes_base(
    rich_for_slope,
    slope_filter_row
  )
  
  list(
    rich_for_slope = rich_for_slope,
    slope = bind_rows(perc, allc, bal$slope),
    raw_resamples = if (save_raw_resamples) bal$raw else NULL
  )
}
run_one_qc_label <- function(qc_row, rich_for_slope_by_filter) {
  
  qc_row <- as_tibble(qc_row)
  
  filter_id <- qc_row$slope_filter_id[[1]]
  qc_name_use <- qc_row$qc_name[[1]]
  k_use <- qc_row$min_cells_per_latbin[[1]]
  
  if (!filter_id %in% names(rich_for_slope_by_filter)) {
    stop(
      "Cannot find slope_filter_id = ", filter_id,
      " in rich_for_slope_by_filter. Available filters are: ",
      paste(names(rich_for_slope_by_filter), collapse = ", ")
    )
  }
  
  message("Calculating QC label for: ", qc_name_use)
  
  rich_for_slope <- rich_for_slope_by_filter[[filter_id]]
  
  rich_for_qc <- filter_valid_latbins_for_qc(
    rich_for_slope,
    k_use
  )
  
  hemi_qc <- classify_stage_hemisphere_qc(rich_for_qc, qc_row,rich_for_slope)
  hemi_qc_full <- complete_hemi_qc(hemi_qc, qc_row)
  
  list(
    rich_for_qc = rich_for_qc %>%
      mutate(
        qc_name = qc_name_use,
        slope_filter_id = filter_id,
        data_role = "QC"
      ),
    hemi_qc = hemi_qc_full
  )
}

slope_outputs <- slope_filters %>%
  split(seq_len(nrow(.))) %>%
  map(run_one_slope_filter)

rich_for_slope_QC <- map(slope_outputs, "rich_for_slope") %>%
  bind_rows()

LDG_slope_by_filter <- map(slope_outputs, "slope") %>%
  bind_rows()

rich_for_slope_by_filter <- split(
  rich_for_slope_QC,
  rich_for_slope_QC$slope_filter_id
)

message("Available slope filters: ", paste(names(rich_for_slope_by_filter), collapse = ", "))

qc_label_outputs <- qc_scenarios %>%
  split(seq_len(nrow(.))) %>%
  map(
    run_one_qc_label,
    rich_for_slope_by_filter = rich_for_slope_by_filter
  )

rich_for_qc_QC <- map(qc_label_outputs, "rich_for_qc") %>%
  bind_rows()

hemi_QC_summary <- map(qc_label_outputs, "hemi_qc") %>%
  bind_rows()

if (save_raw_resamples) {
  percell_raw_resamples_QC <- map(slope_outputs, "raw_resamples") %>%
    bind_rows()
}

LDG_slope_all_QC <- LDG_slope_by_filter %>%
  left_join(
    hemi_QC_summary,
    by = c("slope_filter_id", "bin_midpoint", "stage", "hemisphere")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label),
    color = ifelse(is.na(color), "Bad hemisphere", color),
    slope_90deg = slope * 90,
    slope_sampled_range = slope * lat_range_diagnostic
  )

message("Rows in LDG_slope_by_filter: ", nrow(LDG_slope_by_filter))
message("Rows in hemi_QC_summary: ", nrow(hemi_QC_summary))
message("Rows in LDG_slope_all_QC: ", nrow(LDG_slope_all_QC))

print(
  LDG_slope_all_QC %>%
    distinct(qc_name, slope_filter_id) %>%
    count(slope_filter_id, name = "n_qc_scenarios")
)

stopifnot(all(c(
  "qc_name",
  "slope_filter_id",
  "label",
  "method_group",
  "slope_metric",
  "slope"
) %in% names(LDG_slope_all_QC)))

write.csv(rich_for_slope_QC, file.path(out_dir, paste0(analysis_tag, "_rich_for_slope_by_filter.csv")), row.names = FALSE)
write.csv(rich_for_qc_QC, file.path(out_dir, paste0(analysis_tag, "_rich_for_QC.csv")), row.names = FALSE)
write.csv(hemi_QC_summary, file.path(out_dir, paste0(analysis_tag, "_hemisphere_QC_summary.csv")), row.names = FALSE)
write.csv(LDG_slope_by_filter, file.path(out_dir, paste0(analysis_tag, "_LDG_slope_by_filter_ONLY.csv")), row.names = FALSE)
write.csv(LDG_slope_all_QC, file.path(out_dir, paste0(analysis_tag, "_LDG_slope_ALL_methods_QC.csv")), row.names = FALSE)

saveRDS(
  LDG_slope_all_QC,
  file.path(out_dir, paste0(analysis_tag, "_LDG_slope_ALL_methods_QC.rds"))
)

if (save_raw_resamples) {
  write.csv(
    percell_raw_resamples_QC,
    file.path(out_dir, paste0(analysis_tag, "_percell_raw_resamples_by_filter.csv")),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------
# 7. Downstream analyses
# -----------------------------------------------------------------------

LDG_slope_all_QC_climate <- LDG_slope_all_QC %>%
  mutate(
    bin_midpoint = safe_numeric(bin_midpoint),
    slope = safe_numeric(slope),
    method_group = as.character(method_group),
    slope_metric = as.character(slope_metric)
  ) %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid")) %>%
  mutate(climate_state = factor(climate_state, levels = climate_levels))

write.csv(
  LDG_slope_all_QC_climate,
  file.path(out_dir, paste0(analysis_tag, "_LDG_slope_ALL_methods_QC_climate.csv")),
  row.names = FALSE
)

climate_pairs <- combn(climate_levels, 2, simplify = FALSE)

wilcox_QC_results <- expand_grid(method_specs, qc_name = qc_scenarios$qc_name) %>%
  pmap_dfr(function(method_group, slope_metric, qc_name) {
    
    df <- LDG_slope_all_QC_climate %>%
      filter(
        method_group == !!method_group,
        slope_metric == !!slope_metric,
        qc_name == !!qc_name,
        label == "good",
        !is.na(slope),
        climate_state %in% climate_levels
      )
    
    map_df(climate_pairs, function(pair) {
      
      g1 <- df %>% filter(climate_state == pair[1]) %>% pull(slope)
      g2 <- df %>% filter(climate_state == pair[2]) %>% pull(slope)
      wt <- safe_wilcox(g1, g2)
      m1 <- safe_median(g1)
      m2 <- safe_median(g2)
      
      tibble(
        method_group = method_group,
        slope_metric = slope_metric,
        qc_name = qc_name,
        group1 = pair[1],
        group2 = pair[2],
        n1 = sum(!is.na(safe_numeric(g1))),
        n2 = sum(!is.na(safe_numeric(g2))),
        median1 = m1,
        median2 = m2,
        median_diff = ifelse(is.na(m1) | is.na(m2), NA_real_, m1 - m2),
        p_value = wt$p_value,
        w_statistic = wt$w_statistic,
        iqr1 = safe_iqr(g1),
        iqr2 = safe_iqr(g2),
        test_note = wt$test_note
      )
    })
  }) %>%
  group_by(method_group, slope_metric, qc_name) %>%
  mutate(
    p_adjusted = if (all(is.na(p_value))) NA_real_ else p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(method_group, slope_metric, qc_name, p_adjusted)

# Run the same climate-state comparisons separately within each hemisphere.
# BH adjustment is performed within each method x QC scenario x hemisphere.
wilcox_QC_results_by_hemisphere <- expand_grid(
  method_specs,
  qc_name = qc_scenarios$qc_name,
  hemisphere = c("Northern", "Southern")
) %>%
  pmap_dfr(function(method_group, slope_metric, qc_name, hemisphere) {
    
    df <- LDG_slope_all_QC_climate %>%
      filter(
        method_group == !!method_group,
        slope_metric == !!slope_metric,
        qc_name == !!qc_name,
        .data$hemisphere == .env$hemisphere,
        label == "good",
        !is.na(slope),
        climate_state %in% climate_levels
      )
    
    map_df(climate_pairs, function(pair) {
      
      g1 <- df %>% filter(climate_state == pair[1]) %>% pull(slope)
      g2 <- df %>% filter(climate_state == pair[2]) %>% pull(slope)
      wt <- safe_wilcox(g1, g2)
      m1 <- safe_median(g1)
      m2 <- safe_median(g2)
      
      tibble(
        method_group = method_group,
        slope_metric = slope_metric,
        qc_name = qc_name,
        hemisphere = hemisphere,
        group1 = pair[1],
        group2 = pair[2],
        n1 = sum(!is.na(safe_numeric(g1))),
        n2 = sum(!is.na(safe_numeric(g2))),
        median1 = m1,
        median2 = m2,
        median_diff = ifelse(is.na(m1) | is.na(m2), NA_real_, m1 - m2),
        p_value = wt$p_value,
        w_statistic = wt$w_statistic,
        iqr1 = safe_iqr(g1),
        iqr2 = safe_iqr(g2),
        test_note = wt$test_note
      )
    })
  }) %>%
  group_by(method_group, slope_metric, qc_name, hemisphere) %>%
  mutate(
    p_adjusted = if (all(is.na(p_value))) NA_real_ else p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(method_group, slope_metric, qc_name, hemisphere, p_adjusted)

retention_summary <- LDG_slope_all_QC %>%
  group_by(
    method_group,
    slope_metric,
    qc_name,
    hemisphere,
    qc_type,
    occurrence_min,
    collection_min,
    min_cells_per_latbin,
    zone_requirement
  ) %>%
  summarise(
    n_total_stage_hemi = n(),
    n_with_slope = sum(!is.na(slope)),
    n_good = sum(label == "good", na.rm = TRUE),
    n_good_with_slope = sum(label == "good" & !is.na(slope), na.rm = TRUE),
    prop_good_with_slope = n_good_with_slope / n_total_stage_hemi,
    .groups = "drop"
  )

climate_n_QC <- LDG_slope_all_QC_climate %>%
  filter(label == "good", !is.na(slope), climate_state %in% climate_levels) %>%
  group_by(method_group, slope_metric, qc_name, climate_state) %>%
  summarise(
    n = n(),
    n_Northern = sum(hemisphere == "Northern"),
    n_Southern = sum(hemisphere == "Southern"),
    median_slope = safe_median(slope),
    iqr_slope = safe_iqr(slope),
    .groups = "drop"
  )
wilcox_QC_results <- wilcox_QC_results %>%
  mutate(
    qc_filter = sub("_k[0-9]+_.*$", "", qc_name),
    qc_filter_order = match(
      qc_filter,
      c("occurrence5", "occurrence10", "collection5", "collection10")
    ),
    k_order = as.integer(sub("^.*_k([0-9]+)_.*$", "\\1", qc_name)),
    zone_order = if_else(
      grepl("_tropical_temperate_polar$", qc_name),
      2L,
      1L
    )
  ) %>%
  arrange(
    method_group,
    slope_metric,
    qc_filter_order,
    k_order,
    zone_order,
    p_adjusted
  ) %>%
  select(-qc_filter, -qc_filter_order, -k_order, -zone_order)

wilcox_QC_results_by_hemisphere <- wilcox_QC_results_by_hemisphere %>%
  mutate(
    hemisphere = factor(hemisphere, levels = c("Northern", "Southern")),
    qc_filter = sub("_k[0-9]+_.*$", "", qc_name),
    qc_filter_order = match(
      qc_filter,
      c("occurrence5", "occurrence10", "collection5", "collection10")
    ),
    k_order = as.integer(sub("^.*_k([0-9]+)_.*$", "\\1", qc_name)),
    zone_order = if_else(
      grepl("_tropical_temperate_polar$", qc_name),
      2L,
      1L
    )
  ) %>%
  arrange(
    method_group,
    slope_metric,
    qc_filter_order,
    k_order,
    zone_order,
    hemisphere,
    p_adjusted
  ) %>%
  select(-qc_filter, -qc_filter_order, -k_order, -zone_order)

write.csv(wilcox_QC_results, file.path(out_dir, paste0(analysis_tag, "_wilcoxon_ALL_methods.csv")), row.names = FALSE)
write.csv(
  wilcox_QC_results_by_hemisphere,
  file.path(out_dir, paste0(analysis_tag, "_wilcoxon_ALL_methods_by_hemisphere.csv")),
  row.names = FALSE
)
write.csv(retention_summary, file.path(out_dir, paste0(analysis_tag, "_retention_summary_ALL_methods.csv")), row.names = FALSE)
write.csv(climate_n_QC, file.path(out_dir, paste0(analysis_tag, "_climate_n_ALL_methods.csv")), row.names = FALSE)

# -----------------------------------------------------------------------
# 8. Plot helpers
# -----------------------------------------------------------------------

make_qc_retention_summary <- function(method_table = method_specs_main) {
  
  LDG_slope_all_QC %>%
    semi_join(method_table, by = c("method_group", "slope_metric")) %>%
    group_by(
      method_group,
      slope_metric,
      qc_name,
      qc_type,
      occurrence_min,
      collection_min,
      min_cells_per_latbin,
      zone_requirement,
      hemisphere
    ) %>%
    summarise(
      n_total_stage_hemi = n(),
      n_with_slope = sum(!is.na(slope)),
      n_good = sum(label == "good", na.rm = TRUE),
      n_good_with_slope = sum(label == "good" & !is.na(slope), na.rm = TRUE),
      prop_good_with_slope = n_good_with_slope / n_total_stage_hemi,
      .groups = "drop"
    ) %>%
    mutate(
      method_label = case_when(
        method_group == "percentile_latbin_OLS" ~ paste0("Percentile lat-bin OLS: ", slope_metric),
        method_group == "per_cell_balanced_resampling_OLS" ~ "Per-cell balanced resampling OLS",
        method_group == "per_cell_all_cells_OLS" ~ "Per-cell all-cell OLS",
        TRUE ~ paste(method_group, slope_metric, sep = " | ")
      ),
      filter_label = case_when(
        qc_type == "occurrence" ~ paste0("occurrence >= ", occurrence_min),
        qc_type == "collection" ~ paste0("collection_no >= ", collection_min),
        TRUE ~ qc_type
      ),
      k_label = paste0("k = ", min_cells_per_latbin),
      zone_label = case_when(
        zone_requirement == "tropical_temperate" ~ "T + Te",
        zone_requirement == "tropical_temperate_polar" ~ "T + Te + P",
        TRUE ~ zone_requirement
      ),
      qc_plot_label = paste0(filter_label, " | ", k_label, " | ", zone_label),
      hemisphere = factor(hemisphere, levels = c("Northern", "Southern")),
      method_label = factor(method_label, levels = unique(method_label))
    )
}

draw_qc_retention_by_hemisphere <- function(method_table = method_specs_main) {
  
  retention_data <- make_qc_retention_summary(method_table)
  
  if (!nrow(retention_data)) {
    message("No QC retention data.")
    return(NULL)
  }
  
  qc_levels <- retention_data %>%
    distinct(
      qc_plot_label,
      qc_type,
      occurrence_min,
      collection_min,
      min_cells_per_latbin,
      zone_requirement
    ) %>%
    mutate(
      qc_type_order = case_when(
        qc_type == "occurrence" ~ 1,
        qc_type == "collection" ~ 2,
        TRUE ~ 99
      ),
      threshold_order = case_when(
        qc_type == "occurrence" ~ occurrence_min,
        qc_type == "collection" ~ collection_min,
        TRUE ~ NA_real_
      ),
      zone_order = case_when(
        zone_requirement == "tropical_temperate" ~ 1,
        zone_requirement == "tropical_temperate_polar" ~ 2,
        TRUE ~ 99
      )
    ) %>%
    arrange(qc_type_order, threshold_order, min_cells_per_latbin, zone_order) %>%
    pull(qc_plot_label)
  
  retention_data <- retention_data %>%
    mutate(qc_plot_label = factor(qc_plot_label, levels = rev(qc_levels)))
  
  write.csv(
    retention_data,
    file.path(out_dir, paste0(analysis_tag, "_QC_retention_by_hemisphere.csv")),
    row.names = FALSE
  )
  
  plot_list <- retention_data %>%
    split(.$method_label) %>%
    imap(function(df, method_lab) {
      
      p <- ggplot(
        df,
        aes(x = n_good_with_slope, y = qc_plot_label, fill = hemisphere)
      ) +
        geom_col(
          position = position_dodge(width = 0.75),
          width = 0.65,
          colour = "black",
          linewidth = 0.25
        ) +
        scale_fill_manual(values = hemi_cols, name = "Hemisphere") +
        labs(
          x = "Number of retained slopes",
          y = "QC scenario",
          title = method_lab
        ) +
        theme_minimal() +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(size = 10, colour = "black"),
          axis.text.y = element_text(size = 8.5, colour = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          plot.margin = margin(10, 20, 10, 10)
        )
      
      print(p)
      
      safe_name <- clean_name(paste0("QC_retention_by_hemisphere_", method_lab))
      
      ggsave(
        file.path(fig_dir, "jpg", paste0(analysis_tag, "_", safe_name, ".jpg")),
        p,
        width = 8,
        height = 8,
        dpi = 300
      )
      
      ggsave(
        file.path(fig_dir, "pdf", paste0(analysis_tag, "_", safe_name, ".pdf")),
        p,
        width = 8,
        height = 8,
        dpi = 300
      )
      
      p
    })
  
  return(plot_list)
}
get_qc_grid_by_retention <- function(method_table = method_specs_main, min_prop_total = 0.5) {
  
  retention_data <- make_qc_retention_summary(method_table)
  
  eligible_qc <- retention_data %>%
    group_by(method_group, slope_metric, qc_name) %>%
    summarise(
      n_total = sum(n_total_stage_hemi, na.rm = TRUE),
      n_retained = sum(n_good_with_slope, na.rm = TRUE),
      prop_retained = n_retained / n_total,
      n_retained_Northern = sum(n_good_with_slope[hemisphere == "Northern"], na.rm = TRUE),
      n_retained_Southern = sum(n_good_with_slope[hemisphere == "Southern"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(prop_retained >= min_prop_total) %>%
    arrange(method_group, slope_metric, desc(prop_retained), qc_name)
  
  write.csv(
    eligible_qc,
    file.path(out_dir, paste0(analysis_tag, "_QC_retention_over_", min_prop_total * 100, "pct.csv")),
    row.names = FALSE
  )
  
  eligible_qc
}

draw_baseline_vs_filter_threshold <- function(method_use, metric_use) {
  
  threshold_qc_names <- qc_scenarios %>%
    filter(min_cells_per_latbin == 1, zone_requirement == "tropical_temperate") %>%
    pull(qc_name)
  
  wide <- LDG_slope_all_QC %>%
    filter(
      method_group == method_use,
      slope_metric == metric_use,
      qc_name %in% threshold_qc_names,
      label == "good"
    ) %>%
    select(qc_name, bin_midpoint, stage, hemisphere, slope) %>%
    pivot_wider(names_from = qc_name, values_from = slope)
  
  if (!baseline_qc %in% names(wide)) {
    message("Baseline not available for ", method_use, " | ", metric_use)
    return(NULL)
  }
  
  compare_cols <- setdiff(intersect(threshold_qc_names, names(wide)), baseline_qc)
  
  if (!length(compare_cols)) {
    message("No threshold-comparison QC columns for ", method_use, " | ", metric_use)
    return(NULL)
  }
  
  df <- wide %>%
    pivot_longer(
      all_of(compare_cols),
      names_to = "qc_name",
      values_to = "slope_compare"
    ) %>%
    rename(slope_baseline = all_of(baseline_qc)) %>%
    filter(!is.na(slope_baseline), !is.na(slope_compare)) %>%
    left_join(
      qc_scenarios %>% select(qc_name, qc_type, occurrence_min, collection_min),
      by = "qc_name"
    ) %>%
    mutate(
      qc_label = case_when(
        qc_type == "occurrence" ~ paste0("occurrence >= ", occurrence_min),
        qc_type == "collection" ~ paste0("collection_no >= ", collection_min),
        TRUE ~ qc_name
      ),
      qc_label = factor(
        qc_label,
        levels = c("occurrence >= 10", "collection_no >= 5", "collection_no >= 10")
      ),
      hemisphere = factor(hemisphere, levels = c("Northern", "Southern"))
    ) %>%
    drop_na(slope_baseline, slope_compare, qc_label, hemisphere)
  
  if (!nrow(df)) {
    message("No paired baseline-threshold data for ", method_use, " | ", metric_use)
    return(NULL)
  }
  
  lab <- df %>%
    group_by(qc_label) %>%
    summarise(
      n = n(),
      rho = safe_cor(slope_baseline, slope_compare),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("rho = ", round(rho, 2), "\n", "n = ", n))
  
  axis_lim <- safe_range_expand(c(df$slope_baseline, df$slope_compare))
  
  p <- ggplot(df, aes(x = slope_baseline, y = slope_compare, fill = hemisphere)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.55, colour = "black") +
    geom_point(shape = 21, size = 1.8, alpha = 0.75, colour = "black", stroke = 0.3) +
    geom_text(
      data = lab,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.15,
      size = 2.8
    ) +
    scale_fill_manual(values = hemi_cols, name = "Hemisphere") +
    facet_wrap(~ qc_label, ncol = 3) +
    coord_equal(xlim = axis_lim, ylim = axis_lim, expand = FALSE) +
    labs(
      x = "Baseline LDG slope",
      y = "Filtered LDG slope",
      title = paste0(method_use, " | ", metric_use)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 10),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 30, 10, 10)
    )
  
  print(p)
  
  safe_name <- clean_name(
    paste(method_use, metric_use, "slope_recalculation_threshold_sensitivity", sep = "_")
  )
  
  ggsave(file.path(fig_dir, "jpg", paste0(analysis_tag, "_", safe_name, ".jpg")), p, width = 10, height = 4.2, dpi = 300)
  ggsave(file.path(fig_dir, "pdf", paste0(analysis_tag, "_", safe_name, ".pdf")), p, width = 10, height = 4.2, dpi = 300)
  
  p
}

draw_boxplot <- function(method_use, metric_use, qc_use = baseline_qc) {
  
  slope_data <- LDG_slope_all_QC_climate %>%
    filter(method_group == method_use, slope_metric == metric_use, qc_name == qc_use) %>%
    mutate(
      slope_type = case_when(
        hemisphere == "Northern" ~ "Northern",
        hemisphere == "Southern" ~ "Southern",
        TRUE ~ NA_character_
      ),
      slope_value = ifelse(label == "bad", NA_real_, slope)
    )
  
  slope_data_filtered <- slope_data %>%
    drop_na(slope_value, climate_state, slope_type) %>%
    filter(climate_state %in% climate_levels) %>%
    mutate(
      climate_state = factor(climate_state, levels = climate_levels),
      slope_type = factor(slope_type, levels = c("Northern", "Southern"))
    )
  
  if (!nrow(slope_data_filtered)) {
    message("No data for boxplot: ", method_use, " | ", metric_use, " | ", qc_use)
    return(NULL)
  }
  
  state_n_labels <- slope_data_filtered %>%
    count(climate_state, name = "count") %>%
    mutate(
      state_with_n = paste0(climate_state, "\n(n=", count, ")"),
      climate_state = factor(climate_state, levels = climate_levels)
    ) %>%
    arrange(climate_state)
  
  slope_data_filtered <- slope_data_filtered %>%
    left_join(state_n_labels %>% select(climate_state, state_with_n), by = "climate_state") %>%
    mutate(state_with_n = factor(state_with_n, levels = state_n_labels$state_with_n))
  
  # slope_flag <- slope_data_filtered %>%
  #   group_by(state_with_n, slope_type) %>%
  #   mutate(
  #     q1 = quantile(slope_value, 0.25, na.rm = TRUE),
  #     q3 = quantile(slope_value, 0.75, na.rm = TRUE),
  #     iqr = q3 - q1,
  #     lower = q1 - 1.5 * iqr,
  #     upper = q3 + 1.5 * iqr,
  #     is_outlier = slope_value < lower | slope_value > upper
  #   ) %>%
  #   ungroup()
  # 
  # outlier_data <- slope_flag %>%
  #   mutate(slope_value = ifelse(is_outlier, slope_value, NA_real_))
  
  y_min_val <- min(slope_data_filtered$slope_value, na.rm = TRUE)
  y_max_val <- max(slope_data_filtered$slope_value, na.rm = TRUE)
  
  if (!is.finite(y_min_val) || !is.finite(y_max_val) || y_min_val == y_max_val) {
    y_min_val <- -1
    y_max_val <- 1
  }
  
  p <- ggplot(
    slope_data_filtered,
    aes(x = state_with_n, y = slope_value, fill = slope_type)
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
      position = position_dodge(width = 0.75)
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
    scale_fill_manual(values = hemi_cols, name = "Hemisphere") +
    coord_cartesian(clip = "off", xlim = c(1, 5), ylim = c(y_min_val, y_max_val)) +
    labs(
      x = "Climate state",
      y = "Slope value",
      title = paste0(method_use, " | ", metric_use, "\n", qc_use)
    ) +
    annotate(
      "text",
      x = 5.81,
      y = (y_max_val + 0.1) / 2,
      label = "Non-modern-type",
      vjust = 0.5,
      hjust = 0.5,
      size = 4.5,
      angle = 270
    ) +
    annotate(
      "text",
      x = 5.81,
      y = (y_min_val - 0.1) / 2,
      label = "Modern-type",
      vjust = 0.5,
      hjust = 0.5,
      size = 4.5,
      angle = 270
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white"),
      plot.margin = margin(10, 30, 10, 10),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  print(p)
  
  safe_name <- clean_name(paste(method_use, metric_use, qc_use, "boxplot", sep = "_"))
  
  ggsave(file.path(fig_dir, "jpg", paste0(analysis_tag, "_", safe_name, ".jpg")), p, width = 7, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "pdf", paste0(analysis_tag, "_", safe_name, ".pdf")), p, width = 7, height = 5, dpi = 300)
  
  p
}

draw_grouped_qc_boxplot_all <- function(method_use, metric_use) {
  
  slope_data <- LDG_slope_all_QC_climate %>%
    filter(
      method_group == method_use,
      slope_metric == metric_use,
      qc_type %in% c("occurrence", "collection"),
      min_cells_per_latbin %in% c(1, 2)
    ) %>%
    mutate(
      slope_type = case_when(
        hemisphere == "Northern" ~ "Northern",
        hemisphere == "Southern" ~ "Southern",
        TRUE ~ NA_character_
      ),
      slope_value = ifelse(label == "bad", NA_real_, slope),
      climate_state = factor(climate_state, levels = climate_levels),
      qc_group = case_when(
        qc_type == "occurrence" ~ paste0("occurrence >= ", occurrence_min),
        qc_type == "collection" ~ paste0("collection_no >= ", collection_min),
        TRUE ~ qc_type
      ),
      qc_group = factor(
        qc_group,
        levels = c(
          "occurrence >= 5",
          "occurrence >= 10",
          "collection_no >= 5",
          "collection_no >= 10"
        )
      ),
      k_label = factor(
        paste0("k = ", min_cells_per_latbin),
        levels = c("k = 1", "k = 2")
      ),
      zone_label = factor(
        zone_requirement,
        levels = c("tropical_temperate", "tropical_temperate_polar"),
        labels = c("tropical + temperate", "tropical + temperate + polar")
      ),
      slope_type = factor(slope_type, levels = c("Northern", "Southern"))
    ) %>%
    drop_na(slope_value, climate_state, slope_type, qc_group, k_label, zone_label)
  
  if (!nrow(slope_data)) {
    message("No data for grouped QC boxplot: ", method_use, " | ", metric_use)
    return(NULL)
  }
  
  sample_counts <- slope_data %>%
    group_by(qc_group, k_label, zone_label, climate_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(n_label = paste0("n=", count))
  
  slope_flag <- slope_data %>%
    group_by(qc_group, k_label, zone_label, climate_state, slope_type) %>%
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
    mutate(slope_value = ifelse(is_outlier, slope_value, NA_real_))
  
  y_min_val <- min(slope_data$slope_value, na.rm = TRUE)
  y_max_val <- max(slope_data$slope_value, na.rm = TRUE)
  
  if (!is.finite(y_min_val) || !is.finite(y_max_val) || y_min_val == y_max_val) {
    y_min_val <- -1
    y_max_val <- 1
  }
  
  y_range <- y_max_val - y_min_val
  y_plot_min <- y_min_val - 0.12 * y_range
  y_plot_max <- y_max_val
  y_n_label <- y_plot_min + 0.04 * y_range
  
  p <- ggplot(slope_flag, aes(x = climate_state, y = slope_value, fill = slope_type)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.1, ymax = 0.1, fill = "lightblue", alpha = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8, linetype = "dashed") +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), linewidth = 0.35) +
    geom_jitter(
      data = subset(slope_flag, !is_outlier),
      shape = 21,
      size = 0.8,
      alpha = 0.55,
      position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.75),
      show.legend = FALSE
    ) +
    geom_point(
      data = outlier_data,
      shape = 23,
      size = 0.8,
      stroke = 0.5,
      color = "black",
      position = position_dodge(width = 0.75),
      show.legend = FALSE
    ) +
    geom_text(
      data = sample_counts,
      aes(x = climate_state, y = y_n_label, label = n_label),
      inherit.aes = FALSE,
      size = 2.5,
      vjust = 0.5
    ) +
    scale_fill_manual(values = hemi_cols, labels = c("Northern", "Southern")) +
    ggh4x::facet_nested(
      rows = vars(qc_group, k_label),
      cols = vars(zone_label),
      scales = "free_y",
      switch = "y"
    ) +
    coord_cartesian(clip = "on", ylim = c(y_plot_min, y_plot_max)) +
    labs(
      x = "Climate state",
      y = "Slope value",
      fill = "Hemisphere",
      title = paste0(method_use, " | ", metric_use)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      ggh4x.facet.nestline = element_blank(),
      axis.title = element_text(size = 14),
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 10),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.direction = "vertical",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = 9.5, face = "bold"),
      strip.placement = "outside",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 30, 10, 10)
    )
  
  print(p)
  
  safe_name <- clean_name(paste(method_use, metric_use, "all_QC_grouped_boxplot", sep = "_"))
  
  ggsave(file.path(fig_dir, "jpg", paste0(analysis_tag, "_", safe_name, ".jpg")), p, width = 7, height = 9, dpi = 300)
  ggsave(file.path(fig_dir, "pdf", paste0(analysis_tag, "_", safe_name, ".pdf")), p, width = 7, height = 9, dpi = 300)
  
  p
}

draw_time_series <- function(method_use, metric_use, qc_use = baseline_qc) {
  
  slope_data <- LDG_slope_all_QC_climate %>%
    filter(method_group == method_use, slope_metric == metric_use, qc_name == qc_use) %>%
    mutate(
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
      slope_lower_plot = ifelse(label == "bad", NA_real_, safe_numeric(slope_lower_95)),
      slope_upper_plot = ifelse(label == "bad", NA_real_, safe_numeric(slope_upper_95))
    )
  
  if (sum(!is.na(slope_data$slope_value)) < 3) {
    message("No time-series data for ", method_use, " | ", metric_use, " | ", qc_use)
    return(NULL)
  }
  
  show_ci <- method_use %in% c(
    "per_cell_balanced_resampling_OLS",
    "per_cell_all_cells_OLS"
    # ,
    # "percentile_latbin_OLS"
  ) &&
    any(!is.na(slope_data$slope_lower_plot)) &&
    any(!is.na(slope_data$slope_upper_plot))
  
  data(periods)
  data(epochs)
  
  major_boundaries <- periods$max_age
  x_max_val <- max(time_bins_use$max_ma, na.rm = TRUE)
  x_min_val <- min(time_bins_use$min_ma, na.rm = TRUE)
  
  if (show_ci) {
    y_min_val <- min(c(slope_data$slope_value, slope_data$slope_lower_plot), na.rm = TRUE) * 1.01
    y_max_val <- max(c(slope_data$slope_value, slope_data$slope_upper_plot), na.rm = TRUE) * 1.2
  } else {
    y_min_val <- min(slope_data$slope_value, na.rm = TRUE) * 1.01
    y_max_val <- max(slope_data$slope_value, na.rm = TRUE) * 1.2
  }
  
  if (!is.finite(y_min_val) || !is.finite(y_max_val) || y_min_val == y_max_val) {
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
      # alpha = 0.6,
      colour = NA
    )
  )
  
  plot_one_hemi <- function(hemi_use, hemi_col, tag_lab, show_x = FALSE) {
    
    hemi_df <- slope_data %>%
      filter(slope_type == hemi_use)
    
    p <- ggplot(hemi_df, aes(x = bin_midpoint, y = slope_value)) +
      climate_shade_layer +
      geom_vline(
        xintercept = major_boundaries,
        color = "black",
        linewidth = 0.4,
        alpha = 0.8
      )
    
    if (!show_x) {
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
      geom_hline(yintercept = 0, color = "black", linewidth = 0.4, linetype = "dashed")
    
    if (show_ci) {
      p <- p +
        # geom_ribbon(
        #   aes(ymin = slope_lower_plot, ymax = slope_upper_plot),
        #   fill = unname(hemi_col),
        #   alpha = 0.22,
        #   colour = NA,
        #   na.rm = TRUE,
        #   show.legend = FALSE
        # )
        geom_linerange(
          aes(ymin = slope_lower_plot, ymax = slope_upper_plot),
          colour = unname(hemi_col), alpha = 0.8, width = 2, linewidth = 0.75, na.rm = TRUE,show.legend = FALSE
        )
    }
    
    p <- p +
      geom_line(linewidth = 1, color = unname(hemi_col), na.rm = TRUE) +
      geom_point(
        shape=21,
        size = 2,
        stroke = 0.55,
        color = "black",
        fill = unname(hemi_col),
        na.rm = TRUE
      ) +
      # scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 21), na.translate = FALSE) +
      geom_rect(
        aes(
          xmin = x_min_val,
          xmax = x_max_val,
          ymin = y_min_val,
          ymax = y_max_val * 1.1
        ),
        color = "black",
        fill = NA,
        linewidth = 1
      ) +
      scale_x_reverse(
        name = ifelse(show_x, "Time (Ma)", ""),
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
        y = "Slope value",
        tag = tag_lab,
        subtitle = paste(hemi_use, "Hemisphere")
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none",
        plot.subtitle = element_text(size = 12, face = "bold"),
        axis.title.x = if (!show_x) element_blank() else element_text(size = 14),
        axis.text.x = if (!show_x) element_blank() else element_text(size = 14)
      )
    
    if (show_x) {
      p <- p +
        coord_geo(xlim =c(x_max_val,0),
                  pos = as.list(rep("bottom", 2)),
                  dat = list("periods", "era"),
                  height = list(unit(1.45, "lines"), unit(1.45, "lines")),
                  #fill="darkgrey",
                  lab_color="black",
                  rot = list(0, 0), 
                  size = list(4, 4), 
                  abbrv =list(TRUE,FALSE))
    }
    
    p
  }
  
  P_north <- plot_one_hemi("Northern", hemi_cols["Northern"], "(a)", FALSE)
  P_south <- plot_one_hemi("Southern", hemi_cols["Southern"], "(b)", TRUE)
  
  climate_bar <- ggplot(
    data.frame(climate_state = factor(climate_levels, levels = climate_levels)),
    aes(x = 0, y = 0, fill = climate_state)
  ) +
    geom_point(shape = 22, size = 5, alpha = 0) +
    scale_fill_manual(values = climate_colors, name = "Climate state", drop = FALSE) +
    guides(
      fill = guide_legend(
        override.aes = list(alpha = 1, shape = 22, size = 5, colour = NA)
      )
    ) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      plot.margin = margin(0, 10, 0, 10)
    )
  
  legend_grob <- cowplot::get_legend(climate_bar)
  
  p_final <- ((P_north / P_south) / legend_grob) +
    plot_layout(heights = c(10, 10, 1.2)) +
    plot_annotation(title = paste0(method_use, " | ", metric_use, "\n", qc_use)) &
    theme(
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.spacing = unit(0, "cm"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )
  
  print(p_final)
  
  safe_name <- clean_name(paste(method_use, metric_use, qc_use, "time_series", sep = "_"))
  
  ggsave(file.path(fig_dir, "time_series_jpg", paste0(analysis_tag, "_", safe_name, ".jpg")), p_final, width = 8, height = 7, dpi = 300, type = "cairo")
  ggsave(file.path(fig_dir, "time_series_pdf", paste0(analysis_tag, "_", safe_name, ".pdf")), p_final, width = 8, height = 7, dpi = 300)
  
  p_final
}

# -----------------------------------------------------------------------
# 9. Main plots
# -----------------------------------------------------------------------

qc_retention_plot <- draw_qc_retention_by_hemisphere(method_specs_main)

slope_recalculation_sensitivity_list <- method_specs_main %>%
  mutate(plot = map2(method_group, slope_metric, draw_baseline_vs_filter_threshold))

grouped_all_QC_boxplot_list <- method_specs_main %>%
  mutate(plot = map2(method_group, slope_metric, draw_grouped_qc_boxplot_all))

qc_to_plot <- c(
  "occurrence5_k1_tropical_temperate",
  "occurrence10_k1_tropical_temperate",
  "collection10_k1_tropical_temperate",
  "occurrence10_k3_tropical_temperate_polar"
)

boxplot_list <- expand_grid(method_specs_main, qc_name = qc_to_plot) %>%
  mutate(plot = pmap(list(method_group, slope_metric, qc_name), draw_boxplot))

time_series_qc_grid <- get_qc_grid_by_retention(
  method_table = method_specs_main,
  min_prop_total = min_time_series_retention
)

time_series_list <- time_series_qc_grid %>%
  mutate(plot = pmap(list(method_group, slope_metric, qc_name), draw_time_series))

# -----------------------------------------------------------------------
# 10. NH-SH diagnostics for retained QC scenarios
# -----------------------------------------------------------------------

nhsh_result_dir <- file.path(out_dir, "NH_SH_diagnostics_MAIN_methods")
nhsh_figure_dir <- file.path(fig_dir, "NH_SH_diagnostics_MAIN_methods")

dir.create(nhsh_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nhsh_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(nhsh_result_dir, "paired"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(nhsh_result_dir, "cor"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(nhsh_figure_dir, "jpg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(nhsh_figure_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)

meta_df <- time_bins %>%
  mutate(
    bin_midpoint = mid_ma,
    stage = interval_name,
    era = case_when(
      bin_midpoint < 66 ~ "Cenozoic",
      bin_midpoint >= 66 & bin_midpoint < 251.902 ~ "Mesozoic",
      TRUE ~ "Paleozoic"
    )
  ) %>%
  select(bin_midpoint, stage, era)

# -----------------------------------------------------------------------
# Sampling-profile dissimilarity by slope filter only
# -----------------------------------------------------------------------

sampling_dissim_filter <- rich_for_slope_QC %>%
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

paired_MAIN_df <- LDG_slope_all_QC %>%
  semi_join(method_specs_main, by = c("method_group", "slope_metric")) %>%
  mutate(
    bin_midpoint = safe_numeric(bin_midpoint),
    slope = safe_numeric(slope),
    k_median = safe_numeric(k_median),
    k_cv = safe_numeric(k_cv),
    k_sample = safe_numeric(k_sample),
    interval_width = safe_numeric(slope_upper_95) - safe_numeric(slope_lower_95)
  ) %>%
  filter(label == "good", hemisphere %in% c("Northern", "Southern")) %>%
  left_join(meta_df, by = c("bin_midpoint", "stage")) %>%
  left_join(
    sampling_dissim_filter,
    by = c("slope_filter_id", "bin_midpoint")
  ) %>%
  select(
    method_group,
    slope_metric,
    qc_name,
    qc_type,
    occurrence_min,
    collection_min,
    min_cells_per_latbin,
    zone_requirement,
    bin_midpoint,
    stage,
    era,
    hemisphere,
    slope,
    k_median,
    k_cv,
    k_sample,
    interval_width,
    sampling_profile_dissim
  ) %>%
  pivot_wider(
    names_from = hemisphere,
    values_from = c(slope, k_median, k_cv, k_sample, interval_width),
    names_sep = "_"
  ) %>%
  filter(!is.na(slope_Northern), !is.na(slope_Southern)) %>%
  mutate(
    slope_diff_abs = abs(slope_Northern - slope_Southern),
    slope_diff_signed = slope_Northern - slope_Southern,
    k_cv_diff = abs(k_cv_Northern - k_cv_Southern),
    k_cv_diff_signed = k_cv_Northern - k_cv_Southern,
    k_median_asym = abs(log((k_median_Northern + 1) / (k_median_Southern + 1))),
    k_median_diff = abs(k_median_Northern - k_median_Southern),
    k_median_diff_signed = k_median_Northern - k_median_Southern,
    k_sample_asym = abs(log((k_sample_Northern + 1) / (k_sample_Southern + 1))),
    k_sample_diff_signed = k_sample_Northern - k_sample_Southern,
    interval_width_asym = abs(log((interval_width_Northern + 1) / (interval_width_Southern + 1))),
    interval_width_diff_signed = interval_width_Northern - interval_width_Southern,
    mean_interval_width = rowMeans(
      cbind(interval_width_Northern, interval_width_Southern),
      na.rm = TRUE
    ),
    era = factor(era, levels = c("Paleozoic", "Mesozoic", "Cenozoic"))
  )

NH_SH_cor_MAIN <- paired_MAIN_df %>%
  group_by(
    method_group,
    slope_metric,
    qc_name,
    qc_type,
    occurrence_min,
    collection_min,
    min_cells_per_latbin,
    zone_requirement
  ) %>%
  group_modify(~ bind_rows(
    cor_one(.x, "slope_Northern", "slope_Southern", "NH slope vs SH slope"),
    cor_one(.x, "sampling_profile_dissim", "slope_diff_abs", "|NH-SH slope| vs sampling-profile dissimilarity"),
    cor_one(.x, "k_cv_diff", "slope_diff_abs", "|NH-SH slope| vs k_cv difference"),
    cor_one(.x, "k_median_asym", "slope_diff_abs", "|NH-SH slope| vs k_median asymmetry"),
    cor_one(.x, "k_median_diff", "slope_diff_abs", "|NH-SH slope| vs k_median difference"),
    cor_one(.x, "sampling_profile_dissim", "slope_diff_signed", "NH-SH slope vs sampling-profile dissimilarity"),
    cor_one(.x, "k_cv_diff_signed", "slope_diff_signed", "NH-SH slope vs signed k_cv difference"),
    cor_one(.x, "k_median_diff_signed", "slope_diff_signed", "NH-SH slope vs signed k_median difference")
  )) %>%
  ungroup()

write.csv(
  sampling_dissim_filter,
  file.path(nhsh_result_dir, "NH_SH_sampling_profile_dissimilarity_by_slope_filter.csv"),
  row.names = FALSE
)
write.csv(paired_MAIN_df, file.path(nhsh_result_dir, "NH_SH_all_paired_MAIN_methods.csv"), row.names = FALSE)
write.csv(NH_SH_cor_MAIN, file.path(nhsh_result_dir, "NH_SH_all_correlations_MAIN_methods.csv"), row.names = FALSE)

nhsh_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.tag = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    aspect.ratio = 1
  )

draw_sampling_dissim_time_series_by_filter <- function(slope_filter_use) {
  
  message("Drawing sampling-profile dissimilarity time series: ", slope_filter_use)
  
  df <- sampling_dissim_filter %>%
    filter(slope_filter_id == slope_filter_use) %>%
    right_join(
      time_bins_use %>%
        transmute(
          bin_midpoint = mid_ma,
          stage = interval_name,
          era = case_when(
            bin_midpoint < 66 ~ "Cenozoic",
            bin_midpoint >= 66 & bin_midpoint < 251.902 ~ "Mesozoic",
            TRUE ~ "Paleozoic"
          )
        ),
      by = "bin_midpoint"
    ) %>%
    mutate(
      slope_filter_id = slope_filter_use,
      era = factor(era, levels = c("Paleozoic", "Mesozoic", "Cenozoic"))
    ) %>%
    arrange(desc(bin_midpoint))
  
  if (sum(!is.na(df$sampling_profile_dissim)) < 3) {
    message("  Skip: fewer than 3 non-NA rows.")
    return(NULL)
  }
  
  data(periods)
  data(epochs)
  
  x_max_val <- max(time_bins_use$max_ma, na.rm = TRUE)
  y_max_val <- max(df$sampling_profile_dissim, na.rm = TRUE)
  
  if (!is.finite(y_max_val) || y_max_val <= 0) {
    y_max_val <- 1
  }
  
  p <- ggplot(
    df,
    aes(x = bin_midpoint, y = sampling_profile_dissim)
  ) +
    geom_vline(
      xintercept = periods$max_age,
      color = "black",
      linewidth = 0.35,
      alpha = 0.75
    ) +
    geom_line(
      linewidth = 0.8,
      color = "black",
      na.rm = FALSE
    ) +
    geom_point(
      aes(fill = era),
      shape = 21,
      size = 2.4,
      stroke = 0.45,
      color = "black",
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = era_cols,
      drop = FALSE,
      name = "Era"
    ) +
    scale_x_reverse(
      limits = c(x_max_val, 0),
      breaks = seq(500, 0, -50),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_max_val),
      breaks = pretty(c(0, y_max_val), n = 5),
      expand = c(0, 0)
    ) +
    coord_geo(xlim =c(x_max_val,0),
              pos = as.list(rep("bottom", 2)),
              dat = list("periods", "era"),
              height = list(unit(1.35, "lines"), unit(1.35, "lines")),
              #fill="darkgrey",
              lab_color="black",
              rot = list(0, 0), 
              # size = list(4, 4), 
              abbrv =list(TRUE,FALSE)) +
    labs(
      x = "Time (Ma)",
      y = "Sampling-profile dissimilarity",
      title = paste0("NH-SH sampling-profile dissimilarity | ", slope_filter_use)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, colour = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box.margin = margin(t = -6, r = 0, b = 0, l = 0),
      legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
      legend.spacing.y = unit(0.05, "cm"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.margin = margin(5, 5, 2, 5)
    )
  
  print(p)
  
  out_jpg <- file.path(
    nhsh_figure_dir,
    "sampling_dissim_time_series_jpg",
    "by_slope_filter"
  )
  
  out_pdf <- file.path(
    nhsh_figure_dir,
    "sampling_dissim_time_series_pdf",
    "by_slope_filter"
  )
  
  dir.create(out_jpg, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_pdf, recursive = TRUE, showWarnings = FALSE)
  
  safe_name <- clean_name(
    paste(
      slope_filter_use,
      "sampling_profile_dissim_time_series",
      sep = "_"
    )
  )
  
  ggsave(
    file.path(out_jpg, paste0(safe_name, ".jpg")),
    p,
    width = 8,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(out_pdf, paste0(safe_name, ".pdf")),
    p,
    width = 8,
    height = 4,
    dpi = 300
  )
  
  p
}

draw_one_nh_sh_main <- function(method_use, metric_use, qc_use) {
  
  message("Drawing NH-SH diagnostics: ", method_use, " / ", metric_use, " / ", qc_use)
  
  paired_df <- paired_MAIN_df %>%
    filter(method_group == method_use, slope_metric == metric_use, qc_name == qc_use)
  
  safe_name <- clean_name(paste(method_use, metric_use, qc_use, sep = "_"))
  
  paired_subdir <- file.path(nhsh_result_dir, "paired", method_use, metric_use)
  cor_subdir <- file.path(nhsh_result_dir, "cor", method_use, metric_use)
  jpg_subdir <- file.path(nhsh_figure_dir, "jpg", method_use, metric_use)
  pdf_subdir <- file.path(nhsh_figure_dir, "pdf", method_use, metric_use)
  
  dir.create(paired_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cor_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(jpg_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdf_subdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    paired_df,
    file.path(paired_subdir, paste0(safe_name, "_paired.csv")),
    row.names = FALSE
  )
  
  if (nrow(paired_df) < 3) {
    message("  Skip: fewer than 3 paired NH-SH rows.")
    return(list(
      summary = tibble(
        method_group = method_use,
        slope_metric = metric_use,
        qc_name = qc_use,
        n_pairs = nrow(paired_df),
        status = "skipped_too_few_pairs"
      ),
      paired_df = paired_df,
      cor_results = tibble(),
      plot = NULL
    ))
  }
  
  cor_results <- bind_rows(
    cor_one(paired_df, "slope_Northern", "slope_Southern", "NH slope vs SH slope"),
    cor_one(paired_df, "sampling_profile_dissim", "slope_diff_abs", "|NH-SH slope| vs sampling-profile dissimilarity"),
    cor_one(paired_df, "k_cv_diff", "slope_diff_abs", "|NH-SH slope| vs k_cv difference"),
    cor_one(paired_df, "k_median_asym", "slope_diff_abs", "|NH-SH slope| vs k_median asymmetry"),
    cor_one(paired_df, "sampling_profile_dissim", "slope_diff_signed", "NH-SH slope vs sampling-profile dissimilarity"),
    cor_one(paired_df, "k_cv_diff_signed", "slope_diff_signed", "NH-SH slope vs signed k_cv difference"),
    cor_one(paired_df, "k_median_diff_signed", "slope_diff_signed", "NH-SH slope vs signed k_median difference"),
     ) %>%
    mutate(
      method_group = method_use,
      slope_metric = metric_use,
      qc_name = qc_use
    ) %>%
    select(method_group, slope_metric, qc_name, everything())
  
  write.csv(
    cor_results,
    file.path(cor_subdir, paste0(safe_name, "_cor.csv")),
    row.names = FALSE
  )
  
  cor_A <- cor_results %>% filter(test == "NH slope vs SH slope")
  cor_B <- cor_results %>% filter(test == "|NH-SH slope| vs sampling-profile dissimilarity")
  
  lab_A <- paste0(
    "\u03c1 = ", fmt_num(cor_A$estimate),
    "\np = ", fmt_num(cor_A$p_value),
    "\nn = ", cor_A$n
  )
  
  lab_B <- paste0(
    "\u03c1 = ", fmt_num(cor_B$estimate),
    "\np = ", fmt_num(cor_B$p_value),
    "\nn = ", cor_B$n
  )
  
  axis_lim <- safe_range_expand(c(paired_df$slope_Northern, paired_df$slope_Southern))
  
  p1 <- ggplot(paired_df, aes(x = slope_Northern, y = slope_Southern, fill = era)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_point(shape = 21, size = 3, stroke = 0.6, colour = "black", alpha = 0.9) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = lab_A,
      hjust = 1.05,
      vjust = -0.2,
      size = 3.5,
      colour = "black"
    ) +
    scale_fill_manual(values = era_cols, drop = FALSE) +
    coord_equal(xlim = axis_lim, ylim = axis_lim, expand = FALSE) +
    labs(
      x = "Northern Hemisphere slope",
      y = "Southern Hemisphere slope",
      fill = "Era",
      tag = "(a)"
    ) +
    nhsh_theme
  
  p2 <- ggplot(paired_df, aes(x = sampling_profile_dissim, y = slope_diff_abs, fill = era)) +
    geom_point(shape = 21, size = 3, stroke = 0.6, colour = "black", alpha = 0.9) +
    scale_fill_manual(values = era_cols, drop = FALSE) +
    labs(
      x = "NH-SH sampling-profile dissimilarity",
      y = "|Northern slope - Southern slope|",
      fill = "Era",
      tag = "(b)"
    ) +
    nhsh_theme
  
  fit_ok <- complete.cases(paired_df$sampling_profile_dissim, paired_df$slope_diff_abs)
  
  if (
    sum(fit_ok) >= 3 &&
    length(unique(paired_df$sampling_profile_dissim[fit_ok])) > 1 &&
    length(unique(paired_df$slope_diff_abs[fit_ok])) > 1
  ) {
    p2 <- p2 +
      geom_smooth(
        aes(group = 1),
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        colour = "black",
        fill = "grey80",
        linewidth = 0.6
      )
  }
  
  p2 <- p2 +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = lab_B,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.5,
      colour = "black"
    )
  
  final_plot <- wrap_plots(p1, p2, ncol = 2, guides = "collect") +
    plot_annotation(title = paste0(method_use, " | ", metric_use, "\n", qc_use)) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  print(final_plot)
  
  ggsave(file.path(jpg_subdir, paste0(safe_name, ".jpg")), final_plot, width = 8, height = 4, dpi = 300)
  ggsave(file.path(pdf_subdir, paste0(safe_name, ".pdf")), final_plot, width = 8, height = 4, dpi = 300)
  
  list(
    summary = tibble(
      method_group = method_use,
      slope_metric = metric_use,
      qc_name = qc_use,
      n_pairs = nrow(paired_df),
      status = "ok"
    ),
    paired_df = paired_df,
    cor_results = cor_results,
    plot = final_plot
  )
}

nhsh_run_grid <- get_qc_grid_by_retention(
  method_table = method_specs_main,
  min_prop_total = min_time_series_retention
)

nhsh_main_results <- pmap(
  list(
    nhsh_run_grid$method_group,
    nhsh_run_grid$slope_metric,
    nhsh_run_grid$qc_name
  ),
  draw_one_nh_sh_main
)

sampling_dissim_time_series_list <- slope_filters %>%
  mutate(
    plot = map(
      slope_filter_id,
      draw_sampling_dissim_time_series_by_filter
    )
  )

nhsh_summary_main <- bind_rows(map(nhsh_main_results, "summary"))
nhsh_paired_main <- bind_rows(map(nhsh_main_results, "paired_df"))
nhsh_cor_main <- bind_rows(map(nhsh_main_results, "cor_results"))

write.csv(
  nhsh_summary_main,
  file.path(nhsh_result_dir, "NH_SH_run_summary_MAIN_methods_retained_QC.csv"),
  row.names = FALSE
)

write.csv(
  nhsh_paired_main,
  file.path(nhsh_result_dir, "NH_SH_all_paired_combined_MAIN_methods_retained_QC.csv"),
  row.names = FALSE
)

write.csv(
  nhsh_cor_main,
  file.path(nhsh_result_dir, "NH_SH_all_correlations_combined_MAIN_methods_retained_QC.csv"),
  row.names = FALSE
)

message("All analyses and selected diagnostic plots finished.")
