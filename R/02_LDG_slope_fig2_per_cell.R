# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_fig2_per_cell.R
# Purpose: Draw per-cell LDG slope figures using precomputed per-cell OLS slopes
# -----------------------------------------------------------------------

# rm(list=ls())

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

out_dir_main<-"./figures"
out_dir_stage<-sprintf(
  "./figures/LDG_slope_combined_per_cell/%s km %squota %s equal_area latitude bins %s",
  params$spacing,params$level,rich_params$n_lat_bins,method_tag
)

dir.create(out_dir_main,recursive=TRUE,showWarnings=FALSE)
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

safe_num<-function(x){
  x<-suppressWarnings(as.numeric(as.character(x)))
  x[is.nan(x)|is.infinite(x)]<-NA_real_
  x
}

# -----------------------------------------------------------------------
# 1.Read data
# -----------------------------------------------------------------------

rich_df<-read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,params$level
))

LDG_slope<-read.csv(slope_file)
raw_slope_file <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS raw resamples.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

LDG_slope_raw <- read.csv(raw_slope_file) %>%
  mutate(
    bin_midpoint = safe_num(bin_midpoint),
    slope = safe_num(slope),
    intercept = safe_num(intercept),
    hemisphere = as.character(hemisphere)
  )
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

stage_lookup<-rich_df%>%
  distinct(bin_midpoint,stage)

# -----------------------------------------------------------------------
# 3.Prepare per-cell fitted lines with uncertainty ribbon
# -----------------------------------------------------------------------

LDG_slope<-LDG_slope%>%
  mutate(
    bin_midpoint=safe_num(bin_midpoint),
    slope=safe_num(slope),
    intercept=safe_num(intercept),
    slope_lower_95=safe_num(slope_lower_95),
    slope_upper_95=safe_num(slope_upper_95),
    hemisphere=as.character(hemisphere),
    label=as.character(label)
  )

if(!"slope_lower_95"%in%names(LDG_slope))LDG_slope$slope_lower_95<-NA_real_
if(!"slope_upper_95"%in%names(LDG_slope))LDG_slope$slope_upper_95<-NA_real_

line_range <- rich_df %>%
  group_by(bin_midpoint, stage, hemisphere) %>%
  summarise(
    x_min = min(abs_lat, na.rm = TRUE),
    x_max = max(abs_lat, na.rm = TRUE),
    .groups = "drop"
  )

qc_lookup <- LDG_slope %>%
  select(bin_midpoint, stage, hemisphere, label) %>%
  distinct()

ols_lines <- LDG_slope_raw %>%
  left_join(
    line_range,
    by = c("bin_midpoint", "stage", "hemisphere")
  ) %>%
  filter(
    is.finite(slope),
    is.finite(intercept),
    is.finite(x_min),
    is.finite(x_max)
  ) %>%
  group_by(bin_midpoint, stage, hemisphere) %>%
  group_modify(~{
    
    x_grid <- seq(
      unique(.x$x_min)[1],
      unique(.x$x_max)[1],
      length.out = 100
    )
    
    # 保持每次重采样的 slope/intercept 配对
    predicted <- outer(
      .x$intercept,
      rep(1, length(x_grid))
    ) + outer(
      .x$slope,
      x_grid
    )
    
    # 中心线使用中位截距和中位斜率，保持直线
    median_intercept <- median(
      .x$intercept,
      na.rm = TRUE
    )
    
    median_slope <- median(
      .x$slope,
      na.rm = TRUE
    )
    
    tibble(
      abs_lat = x_grid,
      
      fitted_values =
        median_intercept +
        median_slope * x_grid,
      
      # 置信范围使用联合重采样预测
      fitted_lower = apply(
        predicted,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      
      fitted_upper = apply(
        predicted,
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),
      
      n_resamples_valid = nrow(.x)
    )
  }) %>%
  ungroup() %>%
  left_join(
    qc_lookup,
    by = c("bin_midpoint", "stage", "hemisphere")
  ) %>%
  mutate(
    label = ifelse(is.na(label), "bad", label),
    hemisphere_mod = ifelse(
      label == "bad",
      "Poor quality",
      hemisphere
    ),
    hemisphere_mod = factor(
      hemisphere_mod,
      levels = c(
        "Northern",
        "Southern",
        "Poor quality"
      )
    )
  )

# -----------------------------------------------------------------------
# 4.Faceted per-cell LDG slope figure
# -----------------------------------------------------------------------

combined_rich_fig<-ggplot(
  rich_df,
  aes(x=abs_lat,y=qD_normalized,color=hemisphere_mod,shape=hemisphere)
)+
  geom_point(alpha=0.65,size=1)+
  geom_ribbon(
    data=filter(ols_lines,label=="good"),
    aes(x=abs_lat,ymin=fitted_lower,ymax=fitted_upper,fill=hemisphere_mod),
    alpha=0.18,
    colour=NA,
    inherit.aes=FALSE
  )+
  geom_line(
    data=ols_lines,
    aes(x=abs_lat,y=fitted_values,color=hemisphere_mod,linetype=hemisphere_mod),
    linewidth=0.9,
    inherit.aes=FALSE
  )+
  scale_color_manual(
    name="Line",
    values=hemi_cols,
    breaks=c("Northern","Southern","Poor quality")
  )+
  scale_fill_manual(
    values=hemi_cols,
    guide="none"
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
    color=guide_legend(
      title="Line",
      override.aes=list(shape=NA,size=2,linetype="solid")
    ),
    shape=guide_legend(
      title="Cell",
      override.aes=list(color="black",size=2)
    ),
    linetype="none"
  )+
  facet_wrap(
    ~reorder(bin_midpoint,-as.numeric(as.character(bin_midpoint))),
    labeller=as_labeller(function(x){
      stage_lookup$stage[match(as.numeric(x),stage_lookup$bin_midpoint)]
    }),
    ncol=6
  )+
  scale_x_continuous(
    limits=c(0,90),
    breaks=c(0,30,60,90),
    expand=c(0,0)
  )+
  scale_y_continuous(
    breaks=function(y){
      max_val<-100
      mid_val<- 50 #ceiling(max_val/2/10)*10
      c(0,mid_val,max_val)
    }
  )+
  coord_cartesian(
    ylim = c(0, 100)
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
    panel.spacing.x=unit(0.5,"lines"),
    panel.spacing.y=unit(0.01,"lines"),
    panel.border=element_rect(color="black",fill=NA,linewidth=1),
    panel.grid=element_blank(),
    axis.ticks=element_line(color="black",linewidth=0.5),
    axis.ticks.length=unit(0.08,"cm"),
    axis.text.x=element_text(size=8,color="black"),
    axis.title.x=element_text(size=12,color="black"),
    axis.text.y=element_text(size=8,color="black"),
    axis.title.y=element_text(size=12,color="black",angle=90),
    plot.title=element_text(size=12,face="bold",hjust=0.5),
    legend.position="bottom",
    legend.box="horizontal"
  )

print(combined_rich_fig)

com_path<-sprintf(
  "./figures/LDG_slope_facet_%s_km_%s_quota_%s_equal_area_latitude_bins_%s.jpg",
  params$spacing,params$level,rich_params$n_lat_bins,gsub(" ","_",method_tag)
)

ggsave(
  com_path,
  combined_rich_fig,
  width=8,
  height=9,
  dpi=300
)

# -----------------------------------------------------------------------
# 5.Per-stage per-cell LDG slope figures
# -----------------------------------------------------------------------

for(stg in unique(rich_df$stage)){

  df_bin<-rich_df%>%
    filter(stage==stg,!is.na(hemisphere_mod))%>%
    mutate(
      hemisphere_mod=factor(
        hemisphere_mod,
        levels=c("Northern","Southern","Poor quality")
      )
    )

  bin<-unique(df_bin$bin_midpoint)

  olsl_data_bin<-ols_lines%>%
    filter(stage==stg,!is.na(hemisphere_mod))%>%
    mutate(
      hemisphere_mod=factor(
        hemisphere_mod,
        levels=c("Northern","Southern","Poor quality")
      )
    )

  color_levels<-unique(as.character(df_bin$hemisphere_mod))
  color_palette<-hemi_cols[color_levels]

  p<-ggplot(
    df_bin,
    aes(x=abs_lat,y=qD_normalized,color=hemisphere_mod,shape=hemisphere)
  )+
    geom_point(alpha=0.7,size=2)+
    geom_ribbon(
      data=filter(olsl_data_bin,label=="good"),
      aes(x=abs_lat,ymin=fitted_lower,ymax=fitted_upper,fill=hemisphere_mod),
      alpha=0.18,
      colour=NA,
      inherit.aes=FALSE
    )+
    geom_line(
      data=olsl_data_bin,
      aes(x=abs_lat,y=fitted_values,color=hemisphere_mod,linetype=hemisphere_mod),
      linewidth=1,
      inherit.aes=FALSE
    )+
    scale_color_manual(
      name="Line",
      values=color_palette,
      drop=FALSE
    )+
    scale_fill_manual(
      values=color_palette,
      guide="none",
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
    scale_y_continuous(
      breaks = c(0, 50, 100)
    ) +
    coord_cartesian(
      ylim = c(0, 100)
    ) +
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
      legend.position="right",
      legend.box="vertical",
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
    width=7,
    height=5,
    dpi=200
  )
}

print(com_path)
print(out_dir_stage)