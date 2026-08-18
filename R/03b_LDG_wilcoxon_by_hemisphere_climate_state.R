# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 03b_LDG_wilcoxon_by_hemisphere_climate_state.R
# Purpose: Perform climate-state pairwise Wilcoxon tests separately for
#          Northern and Southern Hemisphere LDG slopes.
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
})

source("./R/options.R")

# -----------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------

climate_levels <- c(
  "Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse"
)

if (!exists("rich_params") || is.null(rich_params$n_lat_bins)) {
  stop(
    "rich_params$n_lat_bins is not defined. ",
    "Run this script through 000_Main_new.R or define rich_params first."
  )
}

analysis_tag <- sprintf(
  "%skm %squota %s equal-area latitude bins",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

results_dir <- "./results"

input_path <- file.path(
  results_dir,
  paste0(
    analysis_tag,
    " LDG slope per-cell balanced OLS and climate states.csv"
  )
)

output_full_path <- file.path(
  results_dir,
  paste0(
    analysis_tag,
    " per-cell balanced OLS wilcoxon tests by hemisphere.csv"
  )
)

output_3dp_path <- file.path(
  results_dir,
  paste0(
    analysis_tag,
    " per-cell balanced OLS wilcoxon tests by hemisphere 3dp.csv"
  )
)

baseline_method <- "per_cell_balanced_resampling_OLS"
baseline_metric <- "median_resampled_slope"
baseline_qc <- "occurrence5_k1_tropical_temperate"

# -----------------------------------------------------------------------
# 1. Helper functions
# -----------------------------------------------------------------------

safe_num <- function(x) {
  x <- as.numeric(x)
  x[is.nan(x) | is.infinite(x)] <- NA_real_
  x
}

safe_wilcox_less <- function(x, y, min_group_n = 2L) {
  x <- safe_num(x)
  y <- safe_num(y)
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < min_group_n || length(y) < min_group_n) {
    return(list(
      p_value = NA_real_,
      statistic = NA_real_,
      test_note = paste0("Skipped: n < ", min_group_n)
    ))
  }

  wt <- suppressWarnings(
    wilcox.test(x, y, alternative = "less", exact = FALSE)
  )

  list(
    p_value = wt$p.value,
    statistic = as.numeric(wt$statistic),
    test_note = "Tested"
  )
}

format_3dp <- function(x) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = 3))
}

# -----------------------------------------------------------------------
# 2. Read and filter the manuscript baseline data
# -----------------------------------------------------------------------

if (!file.exists(input_path)) {
  stop("Input file does not exist: ", input_path)
}

slope_cli_df <- read.csv(input_path, check.names = FALSE)

required_columns <- c(
  "hemisphere", "slope", "label", "climate_state",
  "method_group", "slope_metric", "qc_name"
)

missing_columns <- setdiff(required_columns, names(slope_cli_df))
if (length(missing_columns) > 0L) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

slope_cli_df_filter <- slope_cli_df %>%
  mutate(
    slope = safe_num(slope),
    hemisphere = trimws(as.character(hemisphere)),
    climate_state = trimws(as.character(climate_state)),
    method_group = trimws(as.character(method_group)),
    slope_metric = trimws(as.character(slope_metric)),
    qc_name = trimws(as.character(qc_name)),
    label = trimws(as.character(label))
  ) %>%
  filter(
    qc_name == baseline_qc,
    method_group == baseline_method,
    slope_metric == baseline_metric,
    climate_state %in% climate_levels,
    hemisphere %in% c("Northern", "Southern"),
    label == "good",
    !is.na(slope)
  ) %>%
  mutate(
    hemisphere = factor(
      hemisphere,
      levels = c("Northern", "Southern")
    ),
    climate_state = factor(climate_state, levels = climate_levels)
  )

if (nrow(slope_cli_df_filter) == 0L) {
  cat(
    "\nNo good observations remained after baseline filtering for ",
    rich_params$n_lat_bins,
    " latitude bins. Empty test rows will be written.\n",
    sep = ""
  )
}

# -----------------------------------------------------------------------
# 3. Pairwise Wilcoxon tests within each hemisphere
# -----------------------------------------------------------------------

climate_pairs <- combn(climate_levels, 2, simplify = FALSE)

wil_results <- map_df(levels(slope_cli_df_filter$hemisphere), function(hemi) {
  hemi_df <- slope_cli_df_filter %>%
    filter(hemisphere == hemi)

  hemi_results <- map_df(climate_pairs, function(pair) {
    g1 <- hemi_df %>%
      filter(climate_state == pair[1]) %>%
      pull(slope) %>%
      safe_num()

    g2 <- hemi_df %>%
      filter(climate_state == pair[2]) %>%
      pull(slope) %>%
      safe_num()

    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]

    n1 <- length(g1)
    n2 <- length(g2)
    med1 <- if (n1 > 0L) median(g1) else NA_real_
    med2 <- if (n2 > 0L) median(g2) else NA_real_
    wt <- safe_wilcox_less(g1, g2, min_group_n = 2L)

    tibble(
      hemisphere = hemi,
      method_group = baseline_method,
      slope_metric = baseline_metric,
      qc_name = baseline_qc,
      alternative = "group1 less than group2",
      group1 = pair[1],
      group2 = pair[2],
      n1 = n1,
      n2 = n2,
      median1 = med1,
      median2 = med2,
      median_diff = med1 - med2,
      p_value = wt$p_value,
      w_statistic = wt$statistic,
      iqr1 = if (n1 > 0L) IQR(g1) else NA_real_,
      iqr2 = if (n2 > 0L) IQR(g2) else NA_real_,
      test_note = wt$test_note
    )
  }) %>%
    # Match the supplied workflow: BH adjustment of the pairwise tests.
    # Adjustment is performed separately within each hemisphere.
    mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adjusted, p_value)

  hemi_results
}) %>%
  mutate(
    hemisphere = factor(
      hemisphere,
      levels = c("Northern", "Southern")
    )
  ) %>%
  arrange(hemisphere, p_adjusted, p_value)

# -----------------------------------------------------------------------
# 4. Write full-precision and fixed-three-decimal CSV files
# -----------------------------------------------------------------------

write.csv(wil_results, output_full_path, row.names = FALSE, na = "")

decimal_columns <- c(
  "median1", "median2", "median_diff", "p_value",
  "w_statistic", "iqr1", "iqr2", "p_adjusted"
)

wil_results_3dp <- wil_results %>%
  mutate(across(all_of(decimal_columns), format_3dp))

write.csv(
  wil_results_3dp,
  output_3dp_path,
  row.names = FALSE,
  na = "",
  quote = TRUE
)

cat("\nRows retained after baseline filtering:", nrow(slope_cli_df_filter), "\n")
cat("\nCounts by hemisphere and climate state:\n")
print(
  slope_cli_df_filter %>%
    count(hemisphere, climate_state, name = "n") %>%
    complete(hemisphere, climate_state, fill = list(n = 0L))
)
cat("\nWilcoxon results:\n")
print(wil_results, n = Inf)
cat("\nFull-precision CSV:\n", output_full_path, "\n", sep = "")
cat("\nThree-decimal CSV:\n", output_3dp_path, "\n", sep = "")
