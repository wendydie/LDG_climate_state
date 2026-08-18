# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_fig_per_cell.R
# Purpose: Draw per-cell LDG slope figures using precomputed per-cell OLS slopes
# -----------------------------------------------------------------------

# rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(palaeoverse)
  library(cowplot)
  library(patchwork)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

# -----------------------------------------------------------------------
# 0.Settings
# -----------------------------------------------------------------------
set.seed(123)
occurrence_min<-5

# choose: "balanced" or "all_cells"
slope_method<-"balanced"

method_tag<-case_when(
  slope_method=="balanced"~"per-cell balanced OLS",
  slope_method=="all_cells"~"per-cell all-cells OLS",
  TRUE~slope_method
)

slope_file<-case_when(
  slope_method=="balanced"~sprintf(
    "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS.csv",
    params$spacing,params$level,rich_params$n_lat_bins
  ),
  slope_method=="all_cells"~sprintf(
    "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell all-cells OLS.csv",
    params$spacing,params$level,rich_params$n_lat_bins
  )
)

out_dir_all<-"./figures"
out_dir_stage<-sprintf(
  "./figures/LDG slope per stage/%s km %squota %s equal_area latitude bins %s",
  params$spacing,params$level,rich_params$n_lat_bins,method_tag
)

dir.create(out_dir_all,recursive=TRUE,showWarnings=FALSE)
dir.create(out_dir_stage,recursive=TRUE,showWarnings=FALSE)

hemi_cols<-c(
  "Northern"="#0072B2",
  "Southern"="#E69F00",
  "Poor quality"="#D3D3D3"
)

hemi_shapes<-c(
  "Northern"=16,
  "Southern"=17
)

# -----------------------------------------------------------------------
# 1.Read data
# -----------------------------------------------------------------------

rich_df<-read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,params$level
))

LDG_slope<-read.csv(slope_file)

time_bins<-readRDS("./data/time_bins.RDS")

lat_bins<-palaeoverse::lat_bins_area(n=rich_params$n_lat_bins)%>%
  arrange(min)

lat_zone_lookup <- lat_bins %>%
  mutate(
    lat_bin_mid = mid,
    abs_lat_bin_mid = round(abs(mid), 6),
    lat_zone = case_when(
      abs_lat_bin_mid < 30 ~ "tropical",
      abs_lat_bin_mid < 60 ~ "temperate",
      abs_lat_bin_mid <= 90 ~ "polar",
      TRUE ~ NA_character_
    )
  ) %>%
  select(lat_bin_mid, bin, abs_lat_bin_mid, lat_zone)
# -----------------------------------------------------------------------
# 2.Prepare cell-level richness data with baseline QC
# -----------------------------------------------------------------------

rich_df<-rich_df%>%
  filter(nT>=occurrence_min,t<=2*nT)%>%
  mutate(stage=time_bins$interval_name[match(bin_midpoint,time_bins$mid_ma)])%>%
  filter(bin_midpoint<=486.8500)%>%
  mutate(
    bin_index=findInterval(cell_lat,vec=c(lat_bins$min,Inf)),
    bin=lat_bins$bin[bin_index],
    abs_lat=abs(cell_lat),
    hemisphere=case_when(
      cell_lat>=0~"Northern",
      cell_lat<0~"Southern",
      TRUE~NA_character_
    ),
    lat_band_mid=floor(abs_lat/30)*30+15
  )%>%
  left_join(lat_zone_lookup,by="bin")%>%
  filter(
    !is.na(bin),
    !is.na(abs_lat),
    !is.na(abs_lat_bin_mid),
    !is.na(lat_zone),
    !is.na(hemisphere)
  )%>%
  group_by(bin_midpoint)%>%
  mutate(qD_normalized=qD*100/max(qD,na.rm=TRUE))%>%
  ungroup()

# -----------------------------------------------------------------------
# Baseline QC: occurrence5_k1_tropical_temperate
# -----------------------------------------------------------------------

adjacent_df <- has_adjacent_bins(rich_df, lat_bins) %>%
  distinct(bin_midpoint, hemisphere, label) %>%
  transmute(
    bin_midpoint,
    hemisphere,
    has_adjacent_tt = label == "good"
  )

qc_label_df <- rich_df %>%
  group_by(bin_midpoint, hemisphere) %>%
  summarise(
    has_tropical = any(lat_zone == "tropical", na.rm = TRUE),
    has_temperate = any(lat_zone == "temperate", na.rm = TRUE),
    n_valid_lat_bins = n_distinct(abs_lat_bin_mid),
    .groups = "drop"
  ) %>%
  left_join(
    adjacent_df,
    by = c("bin_midpoint", "hemisphere")
  ) %>%
  mutate(
    has_adjacent_tt = coalesce(has_adjacent_tt, FALSE),
    label = ifelse(
      has_tropical & has_temperate & has_adjacent_tt,
      "good",
      "bad"
    ),
    hemisphere_mod = ifelse(label == "bad", "Poor quality", hemisphere),
    color = ifelse(label == "bad", "Poor quality", hemisphere)
  )

rich_df <- rich_df %>%
  select(-any_of(c("label", "color", "hemisphere_mod"))) %>%
  left_join(
    qc_label_df %>%
      select(bin_midpoint, hemisphere, label, color, hemisphere_mod),
    by = c("bin_midpoint", "hemisphere")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label),
    color = ifelse(is.na(color), "Poor quality", color),
    hemisphere_mod = ifelse(is.na(hemisphere_mod), "Poor quality", hemisphere_mod)
  )

# -----------------------------------------------------------------------
# 3.Prepare per-cell fitted lines
# -----------------------------------------------------------------------

LDG_slope<-LDG_slope%>%
  mutate(
    bin_midpoint=as.numeric(as.character(bin_midpoint)),
    slope=as.numeric(as.character(slope)),
    intercept=as.numeric(as.character(intercept)),
    hemisphere=as.character(hemisphere),
    label=as.character(label)
  )

LDG_slope <- LDG_slope %>%
  select(-any_of(c("label", "color", "hemisphere_mod"))) %>%
  left_join(
    qc_label_df %>%
      select(bin_midpoint, hemisphere, label, color, hemisphere_mod),
    by = c("bin_midpoint", "hemisphere")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label),
    color = ifelse(is.na(color), "Poor quality", color),
    hemisphere_mod = ifelse(is.na(hemisphere_mod), "Poor quality", hemisphere_mod)
  )

line_range<-rich_df%>%
  group_by(bin_midpoint,stage,hemisphere)%>%
  summarise(
    x_min=min(abs_lat,na.rm=TRUE),
    x_max=max(abs_lat,na.rm=TRUE),
    .groups="drop"
  )

ols_lines<-LDG_slope%>%
  left_join(line_range,by=c("bin_midpoint","stage","hemisphere"))%>%
  mutate(
    hemisphere_mod=ifelse(label=="bad","Poor quality",hemisphere)
  )%>%
  filter(
    !is.na(slope),
    !is.na(intercept),
    !is.na(x_min),
    !is.na(x_max)
  )%>%
  pmap_dfr(function(...){
    z<-tibble(...)
    x_seq<-seq(z$x_min,z$x_max,length.out=100)
    tibble(
      bin_midpoint=z$bin_midpoint,
      stage=z$stage,
      hemisphere=z$hemisphere,
      label=z$label,
      hemisphere_mod=z$hemisphere_mod,
      abs_lat=x_seq,
      fitted_values=z$intercept+z$slope*x_seq
    )
  })

# -----------------------------------------------------------------------
# 4.Faceted per-cell LDG slope figure
# -----------------------------------------------------------------------

LDG_s_plot<-ggplot(
  rich_df,
  aes(x=abs_lat,y=qD_normalized,color=hemisphere_mod,shape=hemisphere)
)+
  geom_point(alpha=0.65,size=1)+
  geom_line(
    data=ols_lines,
    aes(x=abs_lat,y=fitted_values,color=hemisphere_mod,linetype=hemisphere_mod),
    linewidth=0.9,
    inherit.aes=FALSE
  )+
  scale_color_manual(
    name="LDG slope",
    values=hemi_cols,
    breaks=c("Northern","Southern","Poor quality")
  )+
  scale_shape_manual(
    name="Cell",
    values=hemi_shapes,
    breaks=c("Northern","Southern")
  )+
  scale_linetype_manual(
    name="Line",
    values=c(
      "Northern"="solid",
      "Southern"="solid",
      "Poor quality"="solid"
    ),
    breaks=c("Northern","Southern","Poor quality")
  )+
  guides(
    shape=guide_legend(
      title="Cell",
      direction="horizontal",
      nrow=1,
      override.aes=list(color="black",size=2)
    ),
    color=guide_legend(
      title="Line",
      direction="horizontal",
      nrow=1,
      override.aes=list(shape=NA,size=2,linetype="solid")
    ),
    linetype="none"
  )+
  scale_x_continuous(
    limits=c(0,90),
    breaks=c(0,30,60,90),
    expand=c(0,0)
  )+
  scale_y_continuous(
    limits=c(0,100),
    breaks=function(y){
      max_val<-100
      mid_val<- 50 #ceiling(max_val/2/10)*10
      c(0,mid_val,max_val)
    }
  )+
  facet_wrap(
    ~reorder(bin_midpoint,-as.numeric(as.character(bin_midpoint))),
    labeller=as_labeller(function(x){
      paste0(rich_df$stage[match(x,rich_df$bin_midpoint)])
    }),
    ncol=6
  )+
  labs(
    x="Absolute palaeolatitude (°)",
    y="Normalized generic richness",
    title=method_tag
  )+
  theme_minimal()+
  theme(
    strip.text=element_text(size=8,face="bold",margin=margin(1,1,1,1)),
    strip.placement="inside",
    panel.spacing=unit(0.01,"lines"),
    panel.spacing.x=unit(0.5,"lines"),
    panel.border=element_rect(color="black",fill=NA,linewidth=1),
    panel.grid=element_blank(),
    axis.ticks=element_line(color="black",linewidth=0.6),
    axis.ticks.length=unit(0.08,"cm"),
    axis.text.x=element_text(size=8,color="black"),
    axis.title.x=element_text(size=12,color="black"),
    axis.text.y=element_text(size=8,color="black"),
    axis.title.y=element_text(size=12,color="black",angle=90),
    plot.title=element_text(size=12,face="bold",hjust=0.5),
    legend.position="bottom",
    legend.box="vertical",
    legend.spacing.y=unit(-0.05,"pt")
  )

# -----------------------------------------------------------------------
# Geological era labels for facet rows
# -----------------------------------------------------------------------

ncol_facet <- 6

facet_order <- rich_df %>%
  distinct(bin_midpoint, stage) %>%
  mutate(
    bin_midpoint = as.numeric(as.character(bin_midpoint))
  ) %>%
  arrange(desc(bin_midpoint)) %>%
  mutate(
    facet_id = row_number(),
    facet_row = ceiling(facet_id / ncol_facet),
    
    era = case_when(
      bin_midpoint > 251.902 ~ "Palaeozoic",
      bin_midpoint > 66.0    ~ "Mesozoic",
      TRUE                   ~ "Cenozoic"
    )
  )

n_facet_rows <- max(facet_order$facet_row)

era_position <- facet_order %>%
  group_by(era) %>%
  summarise(
    row_min = min(facet_row),
    row_max = max(facet_row),
    row_mid = (row_min + row_max) / 2,
    .groups = "drop"
  ) %>%
  mutate(
    # convert facet-row position to 0–1 coordinates
    y = 1 - (row_mid - 0.5) / n_facet_rows
  )

p_main<-LDG_s_plot+theme(legend.position="none")

legend_obj<-cowplot::get_legend(
  LDG_s_plot+
    theme(
      legend.position="bottom",
      legend.box.margin=margin(0,0,0,0),
      legend.margin=margin(0,0,0,0),
      legend.title=element_text(size=9),
      legend.text=element_text(size=8),
      legend.key=element_rect(fill=NA,color=NA),
      legend.background=element_rect(fill=NA,color=NA),
      legend.box.background=element_rect(fill=NA,color=NA)
    )
)

p_main_with_legend <- p_main +
  inset_element(
    legend_obj,
    left = 0.63,
    bottom = -0.02,
    right = 0.88,
    top = 0.04,
    clip = FALSE,
    on_top = TRUE
  )
era_panel <- ggplot(
  era_position,
  aes(x = 0.5, y = y, label = era)
) +
  geom_text(
    angle = 270,
    size = 3.5,
    fontface = "bold"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 3, 0, 0)
  )

# LDG_s_plot_final<-p_main+
#   inset_element(
#     legend_obj,
#     left=0.63,
#     bottom=-0.02,
#     right=0.88,
#     top=0.04,
#     clip=FALSE,
#     on_top=TRUE
#   )
LDG_s_plot_final <-
  p_main_with_legend +
  era_panel +
  plot_layout(
    widths = c(1, 0.055)
  )

print(LDG_s_plot_final)

gg_path<-sprintf(
  "./figures/%s km %squota %s equal-area latitude bins LDG slopes figure %s.jpg",
  params$spacing,params$level,rich_params$n_lat_bins,method_tag
)

ggsave(
  gg_path,
  LDG_s_plot_final,
  width=8,
  height=9,
  dpi=300
)

# -----------------------------------------------------------------------
# 5.Per-stage per-cell LDG slope figures
# -----------------------------------------------------------------------

for(stg in unique(rich_df$stage)){

  df_bin<-rich_df%>%filter(stage==stg)
  bin<-unique(df_bin$bin_midpoint)
  olsl_data<-ols_lines%>%filter(stage==stg)

  color_levels<-unique(df_bin$hemisphere_mod)
  color_palette<-hemi_cols[color_levels]

  p<-ggplot(
    df_bin,
    aes(x=abs_lat,y=qD_normalized,color=hemisphere_mod,shape=hemisphere)
  )+
    geom_point(alpha=0.7,size=2)+
    geom_line(
      data=olsl_data,
      aes(x=abs_lat,y=fitted_values,color=hemisphere_mod,linetype=hemisphere_mod),
      linewidth=1,
      inherit.aes=FALSE
    )+
    scale_color_manual(
      name="LDG slope",
      values=color_palette,
      drop=FALSE
    )+
    scale_shape_manual(
      name="Cell",
      values=hemi_shapes,
      breaks=c("Northern","Southern")
    )+
    scale_linetype_manual(
      name="Line",
      values=c(
        "Northern"="solid",
        "Southern"="solid",
        "Poor quality"="solid"
      ),
      drop=FALSE
    )+
    scale_x_continuous(
      limits=c(0,90),
      breaks=c(0,30,60,90),
      expand=c(0,0)
    )+
    guides(
      color=guide_legend(
        override.aes=list(shape=NA,linewidth=1)
      ),
      shape=guide_legend(
        override.aes=list(color="black",size=2)
      ),
      linetype="none"
    )+
    labs(
      title=sprintf("%s LDG slope for %s (%s Ma)",method_tag,stg,bin),
      x="Absolute palaeolatitude (°)",
      y="Normalized generic richness"
    )+
    theme_minimal()+
    theme(
      legend.position="bottom",
      legend.box="horizontal",
      plot.title=element_text(size=14,face="bold",hjust=0.5),
      axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,color="black"),
      panel.grid=element_blank(),
      panel.border=element_rect(color="black",fill=NA,linewidth=1)
    )

  file_name<-sprintf(
    "%s/Richness_vs_absLatitude_Bin_%s_%s.jpg",
    out_dir_stage,
    bin,
    gsub(" ","_",method_tag)
  )

  ggsave(
    file_name,
    p,
    width=6,
    height=5,
    dpi=150
  )
}

print(gg_path)
print(out_dir_stage)
# -----------------------------------------------------------------------
# 6.Faceted LDG percentile curves
# -----------------------------------------------------------------------
stage_lookup <- rich_df %>%
  distinct(bin_midpoint, stage) %>%
  mutate(
    bin_midpoint = as.numeric(as.character(bin_midpoint))
  ) %>%
  arrange(bin_midpoint)
percentiles_use <- c("q50", "q60", "q75", "q90", "q95")

richness_percentiles <- rich_df %>%
  group_by(bin_midpoint, lat_bin_mid) %>%
  summarise(
    q50 = quantile(qD_normalized, 0.50, na.rm = TRUE),
    q60 = quantile(qD_normalized, 0.60, na.rm = TRUE),
    q75 = quantile(qD_normalized, 0.75, na.rm = TRUE),
    q90 = quantile(qD_normalized, 0.90, na.rm = TRUE),
    q95 = quantile(qD_normalized, 0.95, na.rm = TRUE),
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

LDG_fig <- ggplot() +
  geom_point(
    data = rich_df,
    aes(x = cell_lat, y = qD_normalized),
    size = 0.6,
    color = "black",
    alpha = 0.3
  ) +
  geom_line(
    data = richness_percentiles,
    aes(
      x = lat_bin_mid,
      y = richness,
      color = Percentile,
      group = Percentile
    ),
    alpha = 0.7,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~reorder(bin_midpoint, -as.numeric(as.character(bin_midpoint))),
    labeller = as_labeller(function(x) {
      stage_lookup$stage[
        match(as.numeric(x), stage_lookup$bin_midpoint)
      ]
    }),
    ncol = 6
  ) +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = c(-50, 0, 50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Palaeolatitude",
    y = "Normalized generic richness",
    color = "Percentile"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(
      size = 8,
      face = "bold",
      margin = margin(1, 1, 1, 1)
    ),
    strip.placement = "inside",
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.01, "lines"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.position = "bottom"
  )

print(LDG_fig)

LDG_f_path <- sprintf(
  "./figures/LDG_per_stage_facet_%s_km_%s_quota_%s_equal_area_latitude_bins.jpg",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

ggsave(
  LDG_f_path,
  LDG_fig,
  width = 8,
  height = 9,
  dpi = 300
)

print(LDG_f_path)