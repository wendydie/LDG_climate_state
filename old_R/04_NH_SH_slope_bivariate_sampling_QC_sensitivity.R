# -----------------------------------------------------------------------
# Project: LDG_climate_state
# File: 04_NH_SH_slope_bivariate_sampling_QC_sensitivity.R
# Purpose: NH-SH signed slope comparison and sampling diagnostics
#          for ALL QC combinations and ALL quantiles
# Last updated: 2026-05-11
# -----------------------------------------------------------------------

# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
})

source("./R/options.R")

# -----------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------

era_cols <- c(
  "Paleozoic" = "#9BBA7F",
  "Mesozoic"  = "#67C5CA",
  "Cenozoic"  = "#F2D2A2"
)

quantiles_to_plot <- rich_params$percentiles

analysis_tag <- paste0(
  params$spacing, "km_",
  params$level, "quota_",
  rich_params$n_lat_bins, "lat"
)

percentile_tag <- ifelse(
  length(rich_params$percentiles) > 1,
  "allq",
  rich_params$percentiles[1]
)

analysis_tag_qc <- sprintf(
  "%skm_%squota_%slat_%s",
  params$spacing,
  params$level,
  rich_params$n_lat_bins,
  percentile_tag
)

base_result_dir <- file.path(
  "./results/LDG_QC_sensitivity",
  paste0("latbins_", rich_params$n_lat_bins)
)

base_figure_dir <- file.path(
  "./figures/LDG_QC_sensitivity",
  paste0("latbins_", rich_params$n_lat_bins)
)

result_dir <- file.path(
  base_result_dir,
  paste0("NH_SH_diagnostics_", analysis_tag)
)

figure_dir <- file.path(
  base_figure_dir,
  paste0("NH_SH_diagnostics_", analysis_tag)
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(result_dir, "paired"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "cor"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(figure_dir, "jpg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figure_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 2. Read QC outputs
# -----------------------------------------------------------------------

qc_dir <- base_result_dir

LDG_slope_file <- file.path(
  qc_dir,
  paste0(analysis_tag_qc, "_LDG_slope_QC_sensitivity.csv")
)

hemi_QC_file <- file.path(
  qc_dir,
  paste0(analysis_tag_qc, "_hemisphere_QC_summary.csv")
)

latbin_QC_file <- file.path(
  qc_dir,
  paste0(analysis_tag_qc, "_latbin_percentile_summary.csv")
)

if (!file.exists(LDG_slope_file)) {
  stop("Cannot find LDG_slope_QC file:\n", LDG_slope_file)
}

if (!file.exists(hemi_QC_file)) {
  stop("Cannot find hemisphere QC file:\n", hemi_QC_file)
}

if (!file.exists(latbin_QC_file)) {
  stop("Cannot find latbin QC file:\n", latbin_QC_file)
}

LDG_slope_QC <- read.csv(LDG_slope_file)
hemi_QC_summary <- read.csv(hemi_QC_file)
latbin_QC_summary <- read.csv(latbin_QC_file)

time_bins <- readRDS("./data/time_bins.RDS")

need_cols_slope <- c(
  "qc_name",
  "bin_midpoint",
  "slope",
  "quantile",
  "hemisphere",
  "label",
  "k_median",
  "k_cv"
)

miss_slope <- setdiff(need_cols_slope, names(LDG_slope_QC))

if (length(miss_slope) > 0) {
  stop("Missing columns in LDG_slope_QC: ", paste(miss_slope, collapse = ", "))
}

need_cols_latbin <- c(
  "qc_name",
  "bin_midpoint",
  "hemisphere",
  "abs_lat_bin_mid",
  "n_cells_latbin"
)

miss_latbin <- setdiff(need_cols_latbin, names(latbin_QC_summary))

if (length(miss_latbin) > 0) {
  stop("Missing columns in latbin_QC_summary: ", paste(miss_latbin, collapse = ", "))
}

# -----------------------------------------------------------------------
# 3. Metadata
# -----------------------------------------------------------------------

meta_df <- time_bins %>%
  mutate(
    bin_midpoint = mid_ma,
    stage = interval_name,
    era = case_when(
      bin_midpoint < 66 ~ "Cenozoic",
      bin_midpoint < 251.902 ~ "Mesozoic",
      TRUE ~ "Paleozoic"
    )
  ) %>%
  select(bin_midpoint, stage, era)

# -----------------------------------------------------------------------
# 4. Sampling-profile dissimilarity
# -----------------------------------------------------------------------

sampling_profile_df <- latbin_QC_summary %>%
  transmute(
    qc_name = as.character(qc_name),
    bin_midpoint = as.numeric(bin_midpoint),
    hemisphere = as.character(hemisphere),
    abs_lat_bin_mid = as.numeric(abs_lat_bin_mid),
    n_cells_latbin = as.numeric(n_cells_latbin)
  ) %>%
  distinct(
    qc_name,
    bin_midpoint,
    hemisphere,
    abs_lat_bin_mid,
    .keep_all = TRUE
  ) %>%
  pivot_wider(
    names_from = hemisphere,
    values_from = n_cells_latbin,
    values_fill = 0
  ) %>%
  group_by(qc_name, bin_midpoint) %>%
  summarise(
    sampling_profile_dissim = {
      nN <- if ("Northern" %in% names(cur_data())) Northern else numeric(0)
      nS <- if ("Southern" %in% names(cur_data())) Southern else numeric(0)
      
      if (
        length(nN) == 0 ||
        length(nS) == 0 ||
        sum(nN, na.rm = TRUE) <= 0 ||
        sum(nS, na.rm = TRUE) <= 0
      ) {
        NA_real_
      } else {
        pN <- nN / sum(nN, na.rm = TRUE)
        pS <- nS / sum(nS, na.rm = TRUE)
        0.5 * sum(abs(pN - pS), na.rm = TRUE)
      }
    },
    .groups = "drop"
  )

write.csv(
  sampling_profile_df,
  file.path(result_dir, "NH_SH_sampling_profile_dissimilarity.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------
# 5. Helper functions
# -----------------------------------------------------------------------

cor_fun <- function(x, y) {
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 3) {
    return(data.frame(
      n = sum(ok),
      estimate = NA_real_,
      p_value = NA_real_
    ))
  }
  
  z <- suppressWarnings(
    cor.test(
      x[ok],
      y[ok],
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    n = sum(ok),
    estimate = unname(z$estimate),
    p_value = z$p.value
  )
}

cor_one <- function(data, x, y, test_name) {
  cor_fun(data[[x]], data[[y]]) %>%
    mutate(
      test = test_name,
      x_var = x,
      y_var = y
    ) %>%
    select(test, x_var, y_var, n, estimate, p_value)
}

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    format(round(x, digits), nsmall = digits)
  )
}

safe_range_expand <- function(x, expand_frac = 0.08) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(c(-1, 1))
  }
  
  xr <- range(x)
  
  if (xr[1] == xr[2]) {
    xr <- xr + c(-0.5, 0.5)
  } else {
    xr <- xr + c(-1, 1) * diff(xr) * expand_frac
  }
  
  xr
}

base_theme <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
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

# -----------------------------------------------------------------------
# 6. Draw one QC × one quantile
# -----------------------------------------------------------------------

draw_one_qc_nh_sh <- function(qc_use, quantile_use) {
  
  message("Drawing NH-SH diagnostics: ", qc_use, " / ", quantile_use)
  
  paired_df <- LDG_slope_QC %>%
    mutate(
      qc_name = as.character(qc_name),
      bin_midpoint = as.numeric(bin_midpoint),
      slope = as.numeric(slope),
      k_median = as.numeric(k_median),
      k_cv = as.numeric(k_cv),
      quantile = as.character(quantile)
    ) %>%
    filter(
      qc_name == qc_use,
      quantile == quantile_use,
      label == "good",
      hemisphere %in% c("Northern", "Southern")
    ) %>%
    left_join(meta_df, by = "bin_midpoint") %>%
    select(
      qc_name,
      bin_midpoint,
      stage,
      era,
      hemisphere,
      slope,
      k_median,
      k_cv
    ) %>%
    pivot_wider(
      names_from = hemisphere,
      values_from = c(slope, k_median, k_cv),
      names_sep = "_"
    ) %>%
    distinct(qc_name, bin_midpoint, .keep_all = TRUE) %>%
    left_join(
      sampling_profile_df,
      by = c("qc_name", "bin_midpoint")
    ) %>%
    filter(
      !is.na(slope_Northern),
      !is.na(slope_Southern)
    ) %>%
    mutate(
      slope_diff_abs = abs(slope_Northern - slope_Southern),
      slope_diff_signed = slope_Northern - slope_Southern,
      k_cv_diff = abs(k_cv_Northern - k_cv_Southern),
      k_cv_diff_signed = k_cv_Northern - k_cv_Southern,
      k_median_asym = abs(log((k_median_Northern + 1) / (k_median_Southern + 1))),
      k_median_diff_signed = k_median_Northern - k_median_Southern,
      era = factor(era, levels = c("Paleozoic", "Mesozoic", "Cenozoic")),
      quantile = quantile_use
    )
  
  paired_subdir <- file.path(result_dir, "paired", quantile_use)
  cor_subdir <- file.path(result_dir, "cor", quantile_use)
  jpg_subdir <- file.path(figure_dir, "jpg", quantile_use)
  pdf_subdir <- file.path(figure_dir, "pdf", quantile_use)
  
  dir.create(paired_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cor_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(jpg_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdf_subdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    paired_df,
    file.path(paired_subdir, paste0(qc_use, "_paired.csv")),
    row.names = FALSE
  )
  
  if (nrow(paired_df) < 3) {
    message("  Skip: fewer than 3 paired NH-SH rows.")
    
    summary_row <- tibble(
      qc_name = qc_use,
      quantile = quantile_use,
      n_pairs = nrow(paired_df),
      status = "skipped_too_few_pairs"
    )
    
    return(list(
      summary = summary_row,
      paired_df = paired_df,
      cor_results = tibble(),
      plot = NULL
    ))
  }
  
  cor_results <- bind_rows(
    cor_one(
      paired_df,
      "slope_Northern",
      "slope_Southern",
      "NH slope vs SH slope"
    ),
    cor_one(
      paired_df,
      "sampling_profile_dissim",
      "slope_diff_abs",
      "|NH-SH slope| vs sampling-profile dissimilarity"
    ),
    cor_one(
      paired_df,
      "k_cv_diff",
      "slope_diff_abs",
      "|NH-SH slope| vs k_cv difference"
    ),
    cor_one(
      paired_df,
      "k_median_asym",
      "slope_diff_abs",
      "|NH-SH slope| vs k_median asymmetry"
    ),
    cor_one(
      paired_df,
      "sampling_profile_dissim",
      "slope_diff_signed",
      "NH-SH slope vs sampling-profile dissimilarity"
    ),
    cor_one(
      paired_df,
      "k_cv_diff_signed",
      "slope_diff_signed",
      "NH-SH slope vs signed k_cv difference"
    ),
    cor_one(
      paired_df,
      "k_median_diff_signed",
      "slope_diff_signed",
      "NH-SH slope vs signed k_median difference"
    )
  ) %>%
    mutate(
      qc_name = qc_use,
      quantile = quantile_use
    ) %>%
    select(qc_name, quantile, everything())
  
  write.csv(
    cor_results,
    file.path(cor_subdir, paste0(qc_use, "_cor.csv")),
    row.names = FALSE
  )
  
  cor_A <- cor_results %>%
    filter(test == "NH slope vs SH slope")
  
  cor_B <- cor_results %>%
    filter(test == "|NH-SH slope| vs sampling-profile dissimilarity")
  
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
  
  axis_lim <- safe_range_expand(
    c(paired_df$slope_Northern, paired_df$slope_Southern)
  )
  
  p1 <- ggplot(
    paired_df,
    aes(x = slope_Northern, y = slope_Southern, fill = era)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.4
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_point(
      shape = 21,
      size = 3,
      stroke = 0.6,
      colour = "black",
      alpha = 0.9
    ) +
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
    coord_equal(
      xlim = axis_lim,
      ylim = axis_lim,
      expand = FALSE
    ) +
    labs(
      x = "Northern Hemisphere slope",
      y = "Southern Hemisphere slope",
      fill = "Era",
      tag = "A"
    ) +
    base_theme
  
  p2 <- ggplot(
    paired_df,
    aes(
      x = sampling_profile_dissim,
      y = slope_diff_abs,
      fill = era
    )
  ) +
    geom_point(
      shape = 21,
      size = 3,
      stroke = 0.6,
      colour = "black",
      alpha = 0.9
    ) +
    scale_fill_manual(values = era_cols, drop = FALSE) +
    labs(
      x = "NH-SH sampling-profile dissimilarity",
      y = "|Northern slope - Southern slope|",
      fill = "Era",
      tag = "B"
    ) +
    base_theme
  
  fit_ok <- complete.cases(
    paired_df$sampling_profile_dissim,
    paired_df$slope_diff_abs
  )
  
  ok_fit <- sum(fit_ok) >= 3 &&
    length(unique(paired_df$sampling_profile_dissim[fit_ok])) > 1 &&
    length(unique(paired_df$slope_diff_abs[fit_ok])) > 1
  
  if (ok_fit) {
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
  
  final_plot <- wrap_plots(
    p1,
    p2,
    ncol = 2,
    guides = "collect"
  ) +
    plot_annotation(
      title = paste0(qc_use, " (", quantile_use, ")")
    ) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
    )
  
  print(final_plot)
  
  ggsave(
    file.path(jpg_subdir, paste0(qc_use, ".jpg")),
    final_plot,
    width = 8,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(pdf_subdir, paste0(qc_use, ".pdf")),
    final_plot,
    width = 8,
    height = 4,
    dpi = 300
  )
  
  summary_row <- tibble(
    qc_name = qc_use,
    quantile = quantile_use,
    n_pairs = nrow(paired_df),
    status = "ok"
  )
  
  list(
    summary = summary_row,
    paired_df = paired_df,
    cor_results = cor_results,
    plot = final_plot
  )
}

# -----------------------------------------------------------------------
# 7. Run all QC combinations × all quantiles
# -----------------------------------------------------------------------

all_qc <- LDG_slope_QC %>%
  distinct(qc_name) %>%
  pull(qc_name) %>%
  as.character() %>%
  sort()

quantiles_to_plot <- intersect(
  quantiles_to_plot,
  unique(as.character(LDG_slope_QC$quantile))
)

run_grid <- expand.grid(
  qc_name = all_qc,
  quantile = quantiles_to_plot,
  stringsAsFactors = FALSE
)

all_results <- pmap(
  list(run_grid$qc_name, run_grid$quantile),
  draw_one_qc_nh_sh
)

# -----------------------------------------------------------------------
# 8. Save combined summary tables
# -----------------------------------------------------------------------

summary_all <- bind_rows(map(all_results, "summary"))
paired_all <- bind_rows(map(all_results, "paired_df"))
cor_all <- bind_rows(map(all_results, "cor_results"))

write.csv(
  summary_all,
  file.path(result_dir, "NH_SH_run_summary.csv"),
  row.names = FALSE
)

write.csv(
  paired_all,
  file.path(result_dir, "NH_SH_all_paired.csv"),
  row.names = FALSE
)

write.csv(
  cor_all,
  file.path(result_dir, "NH_SH_all_correlations.csv"),
  row.names = FALSE
)

message("All NH-SH QC diagnostic plots finished.")