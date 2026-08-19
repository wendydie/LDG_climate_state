# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 00_data_distribution_map.R
# Last updated: 2025-10-15
# Author: Die (Wendy) Wen
# -----------------------------------------------------------------------
# Load necessary libraries
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)

# Data preparation ------------------------------------------------------
if (params$clean_again || !file.exists("./data/processed/pbdb_data.RDS")){
  source("./R/00_data_preparation.R")
} else {
  occdf <- readRDS("./data/processed/pbdb_data.RDS")
}

time_bins <- readRDS("./data/time_bins.RDS")
# Get the modern world map (landmasses) as an sf object
world <- ne_countries(scale = "medium", returnclass = "sf")
# base_map_path <- sprintf("data/SC16/tmp/static_polygons_reconstructed_%sMa.shp", 250)
# if (!file.exists(base_map_path)) stop(paste("Base map file not found:", base_map_path))
# world  <- st_read(base_map_path)
# Define Mollweide projection
moll_crs <- "+proj=moll"

# Convert world map to Mollweide projection
world_moll <- st_transform(world, crs = moll_crs)

# Process the fossil occurrence dataset, keeping only unique collection_no records
colldf <- unique(occdf[, c("collection_no", "p_lng", "p_lat", "bin_midpoint")])

# Remove missing and invalid values
colldf <- colldf %>%
  filter(!is.na(p_lng) & !is.na(p_lat) & !is.na(bin_midpoint)) %>%
  filter(p_lng >= -180 & p_lng <= 180 & p_lat >= -90 & p_lat <= 90)

# Convert bin_midpoint to numeric if necessary
if (is.factor(colldf$bin_midpoint) || is.character(colldf$bin_midpoint)) {
  colldf$bin_midpoint <- as.numeric(as.character(colldf$bin_midpoint))
}

# Convert colldf to an sf object with WGS84 CRS (EPSG:4326)
colldf_sf <- st_as_sf(colldf, coords = c("p_lng", "p_lat"), crs = 4326)

# Transform fossil points to Mollweide projection
colldf_moll <- st_transform(colldf_sf, crs = moll_crs)

# Plot the map with Mollweide projection and grid lines
coll_map <- ggplot() +
  # Add transformed modern landmasses
  geom_sf(data = world_moll, fill = "gray80", color = "black", size = 0.3) +
  # Add transformed fossil collection points
  geom_sf(data = colldf_moll, aes(color = bin_midpoint), size = 1, alpha = 0.7) +
  # Add graticules (latitude and longitude grid)
  geom_sf(data = st_graticule(crs = moll_crs), color = "gray60", linetype = "dashed", size = 0.2) +
  # Set labels for the plot
  labs(color = "Age (Ma)") +  # Change legend title to reflect Age in Ma
  # Use a continuous color gradient (heatmap style)
  scale_color_viridis_c(option = "inferno", direction = -1, na.value = "grey50") +  
  # Adjust coordinate scaling with Mollweide projection and grid lines
  coord_sf(crs = moll_crs, expand = FALSE) +  # Enables graticule labels
  # Use a minimalistic theme
  theme_minimal() +
  # Reduce spacing and adjust aesthetics
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  # Show x-axis (Longitude)
    axis.title.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_blank(),  # Show y-axis (Latitude)
    axis.title.y = element_text(size = 12, color = "black", angle = 90),  
    legend.position = "bottom",
    plot.margin = margin(10,10,10,10)
  )
print(coll_map)

coll_map_path <- "./figures/Collection distribution map.jpg"
ggsave(coll_map_path, coll_map, width = 10, height = 5, dpi = 300)


# Define the new folder for saving maps----------------------------------
output_folder <- "figures/maps_output"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
# Loop through unique bin_midpoint values
unique_bins <- sort(unique(colldf$bin_midpoint), decreasing = TRUE)  # Sort from newest to oldest

for (bin in unique_bins) {
  bin_int <- as.integer(bin)  # Ensure bin_midpoint is an integer
  if (bin_int == 0){
    base_map_path <- world 
  }else{
    # Load the reconstructed base map from the corresponding shapefile
    base_map_path <- sprintf("data/SC16/tmp/static_polygons_reconstructed_%sMa.shp", bin_int)
    if (!file.exists(base_map_path)) stop(paste("Base map file not found:", base_map_path))
    base_map <- st_read(base_map_path)
  }
  
  # Filter the fossil dataset for the current bin_midpoint
  subset_data <- colldf %>% filter(bin_midpoint == bin)
  if (nrow(subset_data) == 0) next  # Skip if no data points
  # Extract the corresponding stage name from the time_bins dataset
  stage_value <- unique(time_bins$interval_name[round(time_bins$mid_ma, 0) == round(bin, 0)])
  print(stage_value)
  # Convert the filtered dataset to an sf object with WGS84 projection
  subset_sf <- st_as_sf(subset_data, coords = c("p_lng", "p_lat"), crs = 4326)
  base_map_moll <- st_transform(base_map, crs = moll_crs)
  subset_sf_moll <- st_transform(subset_sf, crs = moll_crs)
  # Generate the plot
  map_plot <- ggplot() +
    geom_sf(data = base_map_moll, fill = "gray80", color = "black", size = 0.3) +  # Reconstructed base map
    geom_sf(data = subset_sf_moll, color = "blue", size = 1, alpha = 0.7) +  # Fossil collection points
    # Add graticules (latitude and longitude grid)
    geom_sf(data = st_graticule(crs = moll_crs), color = "gray60", linetype = "dashed", size = 0.2) +
    
    # Set the title with bin_midpoint and stage name
    labs(title = paste0("Fossil Collections at ", bin, " Ma (", stage_value, ")")) +
    # Apply Mollweide projection
    coord_sf(crs = moll_crs, expand = FALSE) +
    # Minimalistic theme and adjustments
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 12, color = "black", angle = 90),
          legend.position = "none",  # Remove legend since all points are blue
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10))
  
  # Define the output filename
  output_file <- file.path(output_folder, paste0("map_", bin, "Ma.jpg"))
  # Save the plot
  ggsave(output_file, plot = map_plot, width = 10, height = 5, dpi = 150)
  print(paste("Saved:", output_file))  # Print progress
}
