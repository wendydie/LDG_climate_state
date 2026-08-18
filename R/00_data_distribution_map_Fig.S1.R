# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 00_data_distribution_map_Fig.S1.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
library(patchwork)
source("./R/options.R")

# Data preparation ------------------------------------------------------
if (params$clean_again || !file.exists("./data/processed/pbdb_data.RDS")) {
  source("./R/00_data_preparation.R")
} else {
  occdf <- readRDS("./data/processed/pbdb_data.RDS")
}

world <- ne_countries(scale = "medium", returnclass = "sf")
moll_crs <- "+proj=moll"
world_moll <- st_transform(world, crs = moll_crs)
graticule_moll <- st_transform(st_graticule(), crs = moll_crs)

# bin_midpoint 转数值
occdf$bin_midpoint <- as.numeric(as.character(occdf$bin_midpoint))

# 提取现代 / 古地理坐标
modern_df <- occdf %>%
  distinct(collection_no, lng, lat, bin_midpoint) %>%
  filter(!is.na(lng), !is.na(lat), !is.na(bin_midpoint),
         between(lng, -180, 180), between(lat, -90, 90))

paleo_df <- occdf %>%
  distinct(collection_no, p_lng, p_lat, bin_midpoint) %>%
  filter(!is.na(p_lng), !is.na(p_lat), !is.na(bin_midpoint),
         between(p_lng, -180, 180), between(p_lat, -90, 90))

# 转 sf 并投影
modern_sf <- st_as_sf(modern_df, coords = c("lng", "lat"), crs = 4326) %>%
  st_transform(crs = moll_crs)

paleo_sf <- st_as_sf(paleo_df, coords = c("p_lng", "p_lat"), crs = 4326) %>%
  st_transform(crs = moll_crs)

# 统一作图函数
make_map <- function(point_sf, panel_label) {
  ggplot() +
    geom_sf(data = world_moll, fill = "gray80", color = "black", linewidth = 0.25) +
    geom_sf(data = graticule_moll, color = "gray60", linetype = "dashed", linewidth = 0.2) +
    geom_sf(data = point_sf, aes(color = bin_midpoint), size = 0.7, alpha = 0.7) +
    scale_color_viridis_c(option = "inferno", direction = -1, name = "Age (Ma)") +
    coord_sf(crs = moll_crs, expand = FALSE) +
    labs(title = panel_label) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      legend.position = "bottom",
      plot.margin = margin(2, 2, 2, 2)
    )
}

# A 和 B
pA <- make_map(modern_sf, "(a)")
pB <- make_map(paleo_sf, "(b)")

# 合并成一个图
coll_map <- pA / pB +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(coll_map)

# 导出
ggsave(
  "./figures/Collection_distribution_map_AB.jpg",
  coll_map,
  width = 7,
  height = 8,
  dpi = 300
)
