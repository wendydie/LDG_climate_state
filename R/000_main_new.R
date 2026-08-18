# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 000_Main_new.R
# Last updated: 2026-5-20
# -----------------------------------------------------------------------

rm(list = ls())
# source("./R/04_LDG_completeness_estimate.R")
# -----------------------------------------------------------------------
# 1. Safe source function
# -----------------------------------------------------------------------

safe_source <- function(file) {
  cat("\n============================================================\n")
  cat("Running:", file, "\n")
  cat("============================================================\n")
  ok <- tryCatch(
    {
      withCallingHandlers(
        {
          source(file)
        },
        warning = function(w) {
          cat("\nWarning in:", file, "\n")
          cat(conditionMessage(w), "\n")
          invokeRestart("muffleWarning")
        }
      )
      
      cat("Finished:", file, "\n")
      TRUE
    },
    error = function(e) {
      cat("\nSkipped due to error in:", file, "\n")
      cat("Error message:\n")
      cat(conditionMessage(e), "\n")
      cat("Continue to next script...\n")
      FALSE
    }
  )
  ok
}
# -----------------------------------------------------------------------
# 2. Settings
# -----------------------------------------------------------------------
percentiles_use <- c("q50", "q60", "q75", "q90", "q95")
lat_bins_to_run <- c(36, 18, 12, 6) #36, 18, 12, 6
script_list <- c(
  # "./R/04_LDG_completeness_estimate_FigS2.R",
  # "./R/02_LDG_slope_per_cell.R",
  # "./R/02_LDG_slope_per_cell_sensitivity_test.R",
  # "./R/02_LDG_slope_fig_per_cell.R",
  # "./R/02_LDG_slope_fig_per_cell_FigS4.R",
  # "./R/02_LDG_percentile_eras_FigS5.R",
  # "./R/02_LDG_slope_fig2_per_cell.R",
  # "./R/02_LDG_slope_fig3_per_cell.R"
  "./R/03_LDG_compared_in_climate_state_per_cell.R",
  # "./R/03b_LDG_wilcoxon_by_hemisphere_climate_state.R"
  "./R/02b_LDG_slope_QC_sensitivity.R",
  # "./R/02_03_LDG_slope_per_cell_allcells_climate_state.R",
  "./R/05_high_latitude_coverage_summary_FigS3.R",
  "./R/06_NH_SH_slope_bivariate_sampling_FigS8.R"
)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 3. Run all latitude-bin settings
# -----------------------------------------------------------------------

run_log_all <- list()

for (n_lat in lat_bins_to_run) {
  
  cat("\n\n############################################################\n")
  cat("Starting analysis for n_lat_bins =", n_lat, "\n")
  cat("############################################################\n\n")
  
  rich_params <- list(
    percentiles = percentiles_use,
    n_lat_bins = n_lat
  )
  
  run_status <- data.frame(
    n_lat_bins = n_lat,
    script = script_list,
    success = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(script_list)) {
    run_status$success[i] <- safe_source(script_list[i])
  }
  
  run_log_all[[as.character(n_lat)]] <- run_status
  
  cat("\n############################################################\n")
  cat("Finished analysis for n_lat_bins =", n_lat, "\n")
  cat("############################################################\n\n")
}

# -----------------------------------------------------------------------
# 4. Save run log
# -----------------------------------------------------------------------

run_log_all <- dplyr::bind_rows(run_log_all)

print(run_log_all)

write.csv(
  run_log_all,
  "./results/000_Main_new_run_status.csv",
  row.names = FALSE
)

cat("\nAll requested latitude-bin analyses finished.\n")