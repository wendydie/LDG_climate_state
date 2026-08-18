# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 01_LDG_calculation.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(dplyr)
library(dggridR)
library(iNEXT)
source("./R/options.R")
source("./R/functions/calculate_Info.R")

# data preparation ------------------------------------------------------
if (params$clean_again || !file.exists("./data/processed/pbdb_data.RDS")){
  source("./R/00_data_preparation.R")
} else {
  occdf <- readRDS("./data/processed/pbdb_data.RDS")
}

# Spatial binning -------------------------------------------------------
# Construct a global grid with cells with approx params$spacing km
dggs <- dgconstruct(spacing = params$spacing, metric=TRUE, resround='nearest')
# Get cells
occdf$cell <- dgGEO_to_SEQNUM(dggs = dggs, 
                               in_lon_deg = occdf[[params$p_lng]], 
                               in_lat_deg = occdf[[params$p_lat]])$seqnum
# Get coordinates from cells
cellcenter <- dgSEQNUM_to_GEO(dggs, occdf$cell)
occdf$cell_lng <- cellcenter$lon_deg
occdf$cell_lat <- cellcenter$lat_deg

cell_xy <- unique(occdf[c('cell', 'cell_lng', 'cell_lat')])
#iNEXT 
# temporal-spatial binning
occdf$tsbin <- paste0(occdf$bin_midpoint, '_', occdf$cell)
occdf_list <- lapply(unique(occdf$bin_assignment), function(i) subset(occdf, bin_assignment == i))

final_rich_df <- do.call(rbind, lapply(seq_along(occdf_list), function(idx) {
  dat <- occdf_list[[idx]]
  
  process_data(dat,
               q = params$q,
               datatype= params$datatype,
               base = params$base,
               level= params$level,
               nboot = params$nboot,
               req_cols=c("occurrence_no", "collection_no", "reference_no"))
  }))
# Split tsbin into bin_midpoint and cell
split_tsbin <- strsplit(final_rich_df $tsbin, "_")

# Extract bin_midpoint and cell
final_rich_df $bin_midpoint <- sapply(split_tsbin, `[`, 1)  # Extract the first part
final_rich_df $cell <- sapply(split_tsbin, `[`, 2)
final_rich_df  <- merge(final_rich_df , cell_xy, by = "cell", all.x = TRUE)
# Define the directory for saving results
output_dir <- "./results/LDG/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)  # Create the directory if it doesn't exist
}

# Save final_rich_df  as a CSV file
write.csv(final_rich_df , file = paste0(output_dir, params$spacing, "_cell_",params$level,"_richness.csv"), 
          row.names = FALSE) #_nonunique_occurrence
