# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 000_Main.R
# Last updated: 2025-10-15
# Author: Die (Wendy) Wen
# Email: geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
rm(list=ls())
source("./R/04_LDG_completeness_estimate.R")

rich_params <- list(
  # Define the list of percentiles to be calculated
  percentiles = c("q50", "q60", "q75", "q90", "q95"),
  n_lat_bins = 36 # It represents number of equal-area latitude bins. 

)
source("./R/02_LDG_slope.R")
# source("./R/02_LDG_slope_sensitivity_test.R")
source("./R/02_LDG_slope_fig.R")
source("./R/02_LDG_slope_fig2.R")
source("./R/02_LDG_slope_fig3.R")
# source("./R/03_LDG_compared_in_climate_state.R")
# source("./R/03_LDG_sensitivity_test.R")

# source("./R/04_NH_SH_slope_bivariate_sampling_QC_sensitivity.R")

rich_params <- list(
  # Define the list of percentiles to be calculated
  percentiles = c("q50", "q60", "q75", "q90", "q95"),
  n_lat_bins = 18 # It represents number of equal-area latitude bins

)
source("./R/02_LDG_slope.R")
# source("./R/02_LDG_slope_sensitivity_test.R")
source("./R/02_LDG_slope_fig.R")
source("./R/02_LDG_slope_fig2.R")
source("./R/02_LDG_slope_fig3.R")
# source("./R/03_LDG_compared_in_climate_state.R")
# source("./R/03_LDG_sensitivity_test.R")

# source("./R/04_NH_SH_slope_bivariate_sampling_QC_sensitivity.R")

rich_params <- list(
  # Define the list of percentiles to be calculated
  percentiles = c("q50", "q60", "q75", "q90", "q95"),
  n_lat_bins = 12 # It represents number of equal-area latitude bins

)
source("./R/02_LDG_slope.R")
# source("./R/02_LDG_slope_sensitivity_test.R")
source("./R/02_LDG_slope_fig.R")
source("./R/02_LDG_slope_fig2.R")
source("./R/02_LDG_slope_fig3.R")
# source("./R/03_LDG_compared_in_climate_state.R")
# source("./R/03_LDG_sensitivity_test.R")

# source("./R/04_NH_SH_slope_bivariate_sampling_QC_sensitivity.R")

rich_params <- list(
  # Define the list of percentiles to be calculated
  percentiles = c("q50", "q60", "q75", "q90", "q95"),
  n_lat_bins = 6 # It represents number of equal-area latitude bins

)
source("./R/02_LDG_slope.R")
# source("./R/02_LDG_slope_sensitivity_test.R")
source("./R/02_LDG_slope_fig.R")
source("./R/02_LDG_slope_fig2.R")
source("./R/02_LDG_slope_fig3.R")
# source("./R/03_LDG_compared_in_climate_state.R")
# source("./R/03_LDG_sensitivity_test.R")
# 
# source("./R/04_NH_SH_slope_bivariate_sampling_QC_sensitivity.R")
