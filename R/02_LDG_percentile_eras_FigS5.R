# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_percentile_eras_FigS5.R
# Purpose: LDG percentile curves separated by geological eras
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(palaeoverse)
  library(grid)
})

source("./R/options.R")

# -----------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------

rich_df <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,
  params$level
))

time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(
  n = rich_params$n_lat_bins
) %>%
  arrange(min)

rich_df <- rich_df %>%
  filter(
    nT >= 5,
    t <= 2 * nT
  ) %>%
  mutate(
    stage = time_bins$interval_name[
      match(bin_midpoint, time_bins$mid_ma)
    ],
    bin_index = findInterval(
      cell_lat,
      c(lat_bins$min, Inf)
    ),
    bin = lat_bins$bin[bin_index]
  ) %>%
  filter(
    bin_midpoint <= 486.85,
    !is.na(stage),
    !is.na(bin)
  ) %>%
  left_join(
    lat_bins %>%
      select(
        bin,
        lat_bin_mid = mid
      ),
    by = "bin"
  ) %>%
  group_by(bin_midpoint) %>%
  mutate(
    qD_normalized = qD * 100 / max(qD, na.rm = TRUE)
  ) %>%
  ungroup()


# -----------------------------------------------------------------------
# 2. Era assignment
# -----------------------------------------------------------------------

stage_order <- rich_df %>%
  distinct(
    bin_midpoint,
    stage
  ) %>%
  arrange(desc(bin_midpoint)) %>%
  mutate(
    id = row_number()
  )

req <- c(
  "Tremadocian",
  "Changhsingian",
  "Induan",
  "Maastrichtian",
  "Danian"
)

miss <- setdiff(
  req,
  stage_order$stage
)

if (length(miss) > 0) {
  stop(
    "Missing stage(s): ",
    paste(miss, collapse = ", ")
  )
}

i1 <- match(
  "Tremadocian",
  stage_order$stage
)

i2 <- match(
  "Changhsingian",
  stage_order$stage
)

i3 <- match(
  "Induan",
  stage_order$stage
)

i4 <- match(
  "Maastrichtian",
  stage_order$stage
)

i5 <- match(
  "Danian",
  stage_order$stage
)

stage_order <- stage_order %>%
  mutate(
    era = case_when(
      between(id, i1, i2) ~ "Palaeozoic",
      between(id, i3, i4) ~ "Mesozoic",
      id >= i5 ~ "Cenozoic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(era)
  )

rich_df <- rich_df %>%
  left_join(
    stage_order %>%
      select(
        bin_midpoint,
        era
      ),
    by = "bin_midpoint"
  ) %>%
  filter(
    !is.na(era)
  )


# -----------------------------------------------------------------------
# 3. Percentiles
# -----------------------------------------------------------------------

percentiles_use <- c(
  "q50",
  "q60",
  "q75",
  "q90",
  "q95"
)

richness_percentiles <- rich_df %>%
  group_by(
    bin_midpoint,
    stage,
    era,
    lat_bin_mid
  ) %>%
  summarise(
    q50 = quantile(
      qD_normalized,
      0.50,
      na.rm = TRUE
    ),
    q60 = quantile(
      qD_normalized,
      0.60,
      na.rm = TRUE
    ),
    q75 = quantile(
      qD_normalized,
      0.75,
      na.rm = TRUE
    ),
    q90 = quantile(
      qD_normalized,
      0.90,
      na.rm = TRUE
    ),
    q95 = quantile(
      qD_normalized,
      0.95,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(percentiles_use),
    names_to = "Percentile",
    values_to = "richness"
  ) %>%
  mutate(
    Percentile = factor(
      Percentile,
      levels = percentiles_use
    )
  )


# -----------------------------------------------------------------------
# 4. Plot function
# -----------------------------------------------------------------------

make_era_plot <- function(
    era_name,
    show_x = FALSE
) {
  
  lev <- stage_order %>%
    filter(
      era == era_name
    ) %>%
    arrange(
      desc(bin_midpoint)
    ) %>%
    pull(stage)
  
  pts <- rich_df %>%
    filter(
      era == era_name
    ) %>%
    mutate(
      stage = factor(
        stage,
        levels = lev
      )
    )
  
  ln <- richness_percentiles %>%
    filter(
      era == era_name
    ) %>%
    mutate(
      stage = factor(
        stage,
        levels = lev
      )
    )
  
  ggplot() +
    
    geom_point(
      data = pts,
      aes(
        x = cell_lat,
        y = qD_normalized
      ),
      size = 0.45,
      color = "black",
      alpha = 0.25
    ) +
    
    geom_line(
      data = ln,
      aes(
        x = lat_bin_mid,
        y = richness,
        color = Percentile,
        group = Percentile
      ),
      linewidth = 0.7,
      alpha = 0.8
    ) +
    
    facet_wrap(
      ~stage,
      ncol = 6,
      drop = FALSE
    ) +
    
    scale_x_continuous(
      breaks = c(
        -50,
        0,
        50
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    coord_cartesian(
      xlim = c(
        -90,
        90
      ),
      ylim = c(
        0,
        100
      )
    ) +
    
    scale_y_continuous(
      breaks = c(
        0,
        50,
        100
      ),
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    
    labs(
      x = if (
        show_x
      ) {
        "Palaeolatitude (°)"
      } else {
        NULL
      },
      y = NULL,
      color = "Percentile"
    ) +
    
    theme_minimal() +
    
    theme(
      
      strip.text = element_text(
        size = 8,
        face = "bold",
        margin = margin(
          1,
          1,
          1,
          1
        )
      ),
      
      strip.placement = "inside",
      
      panel.spacing.x = unit(
        0.5,
        "lines"
      ),
      
      panel.spacing.y = unit(
        0.08,
        "lines"
      ),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      panel.grid = element_blank(),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.5
      ),
      
      axis.ticks.length = unit(
        0.06,
        "cm"
      ),
      
      axis.text = element_text(
        size = 8,
        color = "black"
      ),
      
      axis.title.x = element_text(
        size = 12,
        color = "black"
      ),
      
      legend.position = "none",
      
      # dashed boundary around each era block
      plot.background = element_rect(
        color = "grey40",
        fill = NA,
        linewidth = 0.9,
        linetype = "dashed"
      ),
      
      # Same margin as Fig. S4
      plot.margin = margin(
        7,
        7,
        7,
        7
      )
    )
}


# -----------------------------------------------------------------------
# 5. Era label
# -----------------------------------------------------------------------

era_label <- function(x) {
  
  ggplot() +
    
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = x,
      angle = 270,
      fontface = "bold",
      size = 3.5
    ) +
    
    xlim(
      0,
      1
    ) +
    
    ylim(
      0,
      1
    ) +
    
    theme_void()
}


# -----------------------------------------------------------------------
# 6. Three era blocks
# -----------------------------------------------------------------------

p_pal <- make_era_plot(
  "Palaeozoic",
  FALSE
)

p_mes <- make_era_plot(
  "Mesozoic",
  FALSE
)

p_cen <- make_era_plot(
  "Cenozoic",
  TRUE
)


row_pal <- p_pal +
  era_label(
    "Palaeozoic"
  ) +
  plot_layout(
    widths = c(
      1,
      0.045
    )
  )


row_mes <- p_mes +
  era_label(
    "Mesozoic"
  ) +
  plot_layout(
    widths = c(
      1,
      0.045
    )
  )


row_cen <- p_cen +
  era_label(
    "Cenozoic"
  ) +
  plot_layout(
    widths = c(
      1,
      0.045
    )
  )


# -----------------------------------------------------------------------
# 7. External legend
# -----------------------------------------------------------------------

legend_plot <- ggplot(
  richness_percentiles,
  aes(
    x = lat_bin_mid,
    y = richness,
    color = Percentile,
    group = Percentile
  )
) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      title.position = "left"
    )
  ) +
  
  theme_void() +
  
  theme(
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      size = 9
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.margin = margin(
      0,
      0,
      0,
      0
    ),
    
    legend.box.margin = margin(
      0,
      0,
      0,
      0
    ),
    
    legend.spacing.x = unit(
      2,
      "pt"
    ),
    
    legend.key.height = unit(
      5,
      "pt"
    ),
    
    legend.key.width = unit(
      14,
      "pt"
    )
  )


legend_patch <- wrap_elements(
  full = cowplot::get_legend(
    legend_plot
  )
)


# -----------------------------------------------------------------------
# 8. Combine Fig. S5
# -----------------------------------------------------------------------

gap <- plot_spacer()

n_pal <- ceiling(
  sum(
    stage_order$era == "Palaeozoic"
  ) / 6
)

n_mes <- ceiling(
  sum(
    stage_order$era == "Mesozoic"
  ) / 6
)

n_cen <- ceiling(
  sum(
    stage_order$era == "Cenozoic"
  ) / 6
)


FigS5 <-
  
  row_pal /
  
  gap /
  
  row_mes /
  
  gap /
  
  row_cen /
  
  legend_patch +
  
  plot_layout(
    heights = c(
      n_pal,
      0.12,
      n_mes,
      0.12,
      n_cen,
      0.45
    )
  )


# -----------------------------------------------------------------------
# 9. Common y-axis title
# -----------------------------------------------------------------------

FigS5 <- (
  
  wrap_elements(
    full = textGrob(
      "Normalized generic richness",
      rot = 90,
      gp = gpar(
        fontsize = 12
      )
    )
  ) |
    
    FigS5
  
) +
  
  plot_layout(
    widths = c(
      0.035,
      1
    )
  )


# -----------------------------------------------------------------------
# 10. Save Fig. S5
# -----------------------------------------------------------------------

fig_path <- sprintf(
  "./figures/FigS5_percentile_%s_km_%s_quota_%s_equal-area-latbins.jpg",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

if (
  file.exists(
    fig_path
  )
) {
  file.remove(
    fig_path
  )
}


ggsave(
  filename = fig_path,
  plot = FigS5,
  width = 8,
  height = 9,
  dpi = 300,
  bg = "white"
)


print(
  FigS5
)


cat(
  "\nSaved: ",
  normalizePath(
    fig_path
  ),
  "\n",
  sep = ""
)