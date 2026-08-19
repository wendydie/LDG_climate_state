# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_fig3_per_cell_FigS6S7.R
# Purpose: Fig. S6-S7 hemisphere-separated per-cell LDG slopes by era
# Fig. S6 = Northern Hemisphere
# Fig. S7 = Southern Hemisphere
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(purrr)
  library(palaeoverse); library(cowplot); library(patchwork); library(grid)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

occurrence_min <- 5
slope_method <- "balanced"

pal <- c(
  Northern="#0072B2",
  Southern="#E69F00",
  `Poor quality`="#D3D3D3"
)

hemi_shapes <- c(Northern=16,Southern=17)

safe_num <- function(x){
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[!is.finite(x)] <- NA_real_
  x
}

# -----------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------
rich_df <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,params$level
))

raw_slope_file <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS raw resamples.csv",
  params$spacing,params$level,rich_params$n_lat_bins
)

LDG_slope_raw <- read.csv(raw_slope_file) %>%
  mutate(
    bin_midpoint=safe_num(bin_midpoint),
    slope=safe_num(slope),
    intercept=safe_num(intercept),
    hemisphere=as.character(hemisphere)
  )

time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(
  n=rich_params$n_lat_bins
) %>% arrange(min)

lat_lookup <- lat_bins %>%
  mutate(
    abs_lat_bin_mid=abs(mid),
    lat_zone=case_when(
      abs_lat_bin_mid<30~"tropical",
      abs_lat_bin_mid<60~"temperate",
      abs_lat_bin_mid<=90~"polar",
      TRUE~NA_character_
    )
  ) %>%
  select(bin,abs_lat_bin_mid,lat_zone)

# -----------------------------------------------------------------------
# 2. Richness + QC
# -----------------------------------------------------------------------
rich_df <- rich_df %>%
  filter(nT>=occurrence_min,t<=2*nT) %>%
  mutate(
    stage=time_bins$interval_name[
      match(bin_midpoint,time_bins$mid_ma)
    ]
  ) %>%
  filter(bin_midpoint<=486.85,!is.na(stage)) %>%
  mutate(
    bin_index=findInterval(cell_lat,c(lat_bins$min,Inf)),
    bin=lat_bins$bin[bin_index],
    abs_lat=abs(cell_lat),
    hemisphere=ifelse(cell_lat>=0,"Northern","Southern")
  ) %>%
  left_join(lat_lookup,by="bin") %>%
  filter(
    !is.na(bin),!is.na(abs_lat_bin_mid),
    !is.na(lat_zone),!is.na(hemisphere)
  ) %>%
  group_by(bin_midpoint) %>%
  mutate(qD_normalized=qD*100/max(qD,na.rm=TRUE)) %>%
  ungroup()

adjacent_df <- has_adjacent_bins(rich_df,lat_bins) %>%
  distinct(bin_midpoint,hemisphere,label) %>%
  transmute(
    bin_midpoint,hemisphere,
    has_adjacent_tt=label=="good"
  )

qc <- rich_df %>%
  group_by(bin_midpoint,hemisphere) %>%
  summarise(
    has_tropical=any(lat_zone=="tropical",na.rm=TRUE),
    has_temperate=any(lat_zone=="temperate",na.rm=TRUE),
    .groups="drop"
  ) %>%
  left_join(adjacent_df,by=c("bin_midpoint","hemisphere")) %>%
  mutate(
    has_adjacent_tt=coalesce(has_adjacent_tt,FALSE),
    label=ifelse(
      has_tropical & has_temperate & has_adjacent_tt,
      "good","bad"
    ),
    hemisphere_mod=ifelse(
      label=="bad","Poor quality",hemisphere
    )
  )

rich_df <- rich_df %>%
  left_join(
    qc %>% select(
      bin_midpoint,hemisphere,label,hemisphere_mod
    ),
    by=c("bin_midpoint","hemisphere")
  ) %>%
  mutate(
    label=coalesce(label,"bad"),
    hemisphere_mod=coalesce(
      hemisphere_mod,"Poor quality"
    )
  )

# -----------------------------------------------------------------------
# 3. Geological eras
# -----------------------------------------------------------------------
stage_order <- rich_df %>%
  distinct(bin_midpoint,stage) %>%
  arrange(desc(bin_midpoint)) %>%
  mutate(id=row_number())

required <- c(
  "Tremadocian","Changhsingian",
  "Induan","Maastrichtian","Danian"
)

miss <- setdiff(required,stage_order$stage)
if(length(miss)>0)
  stop("Missing stage(s): ",paste(miss,collapse=", "))

i1 <- match("Tremadocian",stage_order$stage)
i2 <- match("Changhsingian",stage_order$stage)
i3 <- match("Induan",stage_order$stage)
i4 <- match("Maastrichtian",stage_order$stage)
i5 <- match("Danian",stage_order$stage)

stage_order <- stage_order %>%
  mutate(
    era=case_when(
      between(id,i1,i2)~"Palaeozoic",
      between(id,i3,i4)~"Mesozoic",
      id>=i5~"Cenozoic",
      TRUE~NA_character_
    )
  ) %>%
  filter(!is.na(era))

rich_df <- rich_df %>%
  left_join(
    stage_order %>% select(bin_midpoint,era),
    by="bin_midpoint"
  )

# -----------------------------------------------------------------------
# 4. Resampled OLS lines + 95% interval
# -----------------------------------------------------------------------
stage_lookup <- rich_df %>%
  distinct(bin_midpoint,stage)

if(!"stage" %in% names(LDG_slope_raw)){
  LDG_slope_raw <- LDG_slope_raw %>%
    left_join(stage_lookup,by="bin_midpoint")
}

line_range <- rich_df %>%
  group_by(bin_midpoint,stage,hemisphere) %>%
  summarise(
    x_min=min(abs_lat,na.rm=TRUE),
    x_max=max(abs_lat,na.rm=TRUE),
    .groups="drop"
  )

ols_lines <- LDG_slope_raw %>%
  left_join(
    line_range,
    by=c("bin_midpoint","stage","hemisphere")
  ) %>%
  filter(
    is.finite(slope),is.finite(intercept),
    is.finite(x_min),is.finite(x_max)
  ) %>%
  group_by(bin_midpoint,stage,hemisphere) %>%
  group_modify(~{
    
    x <- seq(
      unique(.x$x_min)[1],
      unique(.x$x_max)[1],
      length.out=100
    )
    
    pred <- outer(.x$intercept,rep(1,length(x))) +
      outer(.x$slope,x)
    
    tibble(
      abs_lat=x,
      fitted_values=
        median(.x$intercept,na.rm=TRUE)+
        median(.x$slope,na.rm=TRUE)*x,
      fitted_lower=apply(
        pred,2,quantile,.025,na.rm=TRUE
      ),
      fitted_upper=apply(
        pred,2,quantile,.975,na.rm=TRUE
      )
    )
  }) %>%
  ungroup() %>%
  left_join(
    qc %>% select(
      bin_midpoint,hemisphere,label,hemisphere_mod
    ),
    by=c("bin_midpoint","hemisphere")
  ) %>%
  left_join(
    stage_order %>% select(bin_midpoint,era),
    by="bin_midpoint"
  ) %>%
  mutate(
    label=coalesce(label,"bad"),
    hemisphere_mod=coalesce(
      hemisphere_mod,"Poor quality"
    )
  )

# -----------------------------------------------------------------------
# 5. Era plot
# -----------------------------------------------------------------------
make_era_plot <- function(hemi,era_name,show_x=FALSE){
  
  lev <- stage_order %>%
    filter(era==era_name) %>%
    arrange(desc(bin_midpoint)) %>%
    pull(stage)
  
  df <- rich_df %>%
    filter(
      hemisphere==hemi,
      era==era_name
    ) %>%
    mutate(stage=factor(stage,levels=lev))
  
  ln <- ols_lines %>%
    filter(
      hemisphere==hemi,
      era==era_name
    ) %>%
    mutate(stage=factor(stage,levels=lev))
  
  cols <- setNames(
    c(
      unname(pal[hemi]),
      unname(pal["Poor quality"])
    ),
    c(hemi,"Poor quality")
  )
  
  ggplot(
    df,
    aes(abs_lat,qD_normalized,color=hemisphere_mod)
  )+
    geom_point(alpha=.7,size=.9)+
    geom_ribbon(
      data=filter(ln,label=="good"),
      aes(
        x=abs_lat,
        ymin=fitted_lower,
        ymax=fitted_upper,
        fill=hemisphere_mod
      ),
      alpha=.18,
      colour=NA,
      inherit.aes=FALSE
    )+
    geom_line(
      data=ln,
      aes(
        x=abs_lat,
        y=fitted_values,
        color=hemisphere_mod
      ),
      linewidth=.8,
      inherit.aes=FALSE
    )+
    scale_color_manual(
      values=cols,
      breaks=c(hemi,"Poor quality"),
      drop=FALSE
    )+
    scale_fill_manual(
      values=cols,
      guide="none",
      drop=FALSE
    )+
    scale_x_continuous(
      limits=c(0,90),
      breaks=c(0,30,60,90),
      expand=c(0,0)
    )+
    scale_y_continuous(
      limits=c(0,100),
      breaks=c(0,50,100)
    )+
    facet_wrap(
      ~stage,
      ncol=6,
      drop=FALSE
    )+
    labs(
      x=if(show_x)
        "Absolute palaeolatitude (°)"
      else NULL,
      y=NULL
    )+
    theme_minimal()+
    theme(
      strip.text=element_text(
        size=7,
        face="bold",
        margin=margin(0,0,0,0)
      ),
      strip.placement="inside",
      
      panel.spacing.x=unit(.10,"lines"),
      panel.spacing.y=unit(0,"lines"),
      
      panel.border=element_rect(
        color="black",
        fill=NA,
        linewidth=.7
      ),
      panel.grid=element_blank(),
      
      axis.ticks=element_line(
        color="black",
        linewidth=.4
      ),
      axis.ticks.length=unit(.04,"cm"),
      axis.text=element_text(
        size=7,
        color="black"
      ),
      axis.title.x=element_text(
        size=10,
        color="black",
        margin=margin(t=2)
      ),
      
      legend.position="none",
      
      plot.background=element_rect(
        color="grey40",
        fill=NA,
        linewidth=.75,
        linetype="dashed"
      ),
      
      plot.margin=margin(3,3,3,3)
    )
}

era_label <- function(x){
  ggplot()+
    annotate(
      "text",
      x=.5,y=.5,
      label=x,
      angle=270,
      fontface="bold",
      size=3.3
    )+
    xlim(0,1)+ylim(0,1)+
    theme_void()
}

# -----------------------------------------------------------------------
# 6. One hemisphere figure
# -----------------------------------------------------------------------
make_hemi_fig <- function(hemi,fig_name){
  
  p_pal <- make_era_plot(
    hemi,"Palaeozoic",FALSE
  )
  p_mes <- make_era_plot(
    hemi,"Mesozoic",FALSE
  )
  p_cen <- make_era_plot(
    hemi,"Cenozoic",TRUE
  )
  
  # IMPORTANT:
  # use wrap_plots explicitly to avoid wrap_dims error
  row_pal <- wrap_plots(
    p_pal,
    era_label("Palaeozoic"),
    nrow=1,
    widths=c(1,.025)
  )
  
  row_mes <- wrap_plots(
    p_mes,
    era_label("Mesozoic"),
    nrow=1,
    widths=c(1,.025)
  )
  
  row_cen <- wrap_plots(
    p_cen,
    era_label("Cenozoic"),
    nrow=1,
    widths=c(1,.025)
  )
  
  n_pal <- ceiling(
    sum(stage_order$era=="Palaeozoic")/6
  )
  n_mes <- ceiling(
    sum(stage_order$era=="Mesozoic")/6
  )
  n_cen <- ceiling(
    sum(stage_order$era=="Cenozoic")/6
  )
  
  p_body <- wrap_plots(
    row_pal,
    plot_spacer(),
    row_mes,
    plot_spacer(),
    row_cen,
    ncol=1,
    heights=c(
      n_pal,.015,
      n_mes,.015,
      n_cen
    )
  )
  
  # ---------------------------------------------------------------------
  # Legend outside
  # ---------------------------------------------------------------------
  cols <- setNames(
    c(
      unname(pal[hemi]),
      unname(pal["Poor quality"])
    ),
    c(hemi,"Poor quality")
  )
  
  leg_df <- tibble(
    x=1:2,
    y=1,
    group=factor(
      c(hemi,"Poor quality"),
      levels=c(hemi,"Poor quality")
    )
  )
  
  legend_plot <- ggplot(
    leg_df,
    aes(x,y,color=group)
  )+
    geom_line(
      aes(group=group),
      linewidth=1
    )+
    geom_point(size=1.8)+
    scale_color_manual(
      name="Line",
      values=cols,
      breaks=c(hemi,"Poor quality")
    )+
    guides(
      color=guide_legend(
        nrow=1,
        title.position="left"
      )
    )+
    theme_void()+
    theme(
      legend.position="bottom",
      legend.box="horizontal",
      legend.title=element_text(size=9),
      legend.text=element_text(size=8),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0)
    )
  
  legend_patch <- wrap_elements(
    full=cowplot::get_legend(legend_plot)
  )
  
  p_all <- wrap_plots(
    p_body,
    legend_patch,
    ncol=1,
    heights=c(1,.035)
  )
  
  p_final <- wrap_plots(
    wrap_elements(
      full=textGrob(
        "Normalized generic richness",
        rot=90,
        gp=gpar(fontsize=11)
      )
    ),
    p_all,
    nrow=1,
    widths=c(.022,1)
  )
  
  path <- sprintf(
    "./figures/%s_%skm_%squota_%s_equal-area_latbins_%s.jpg",
    fig_name,
    params$spacing,
    params$level,
    rich_params$n_lat_bins,
    slope_method
  )
  
  if(file.exists(path)) file.remove(path)
  
  ggsave(
    filename=path,
    plot=p_final,
    width=8,
    height=8.5,
    dpi=300,
    bg="white"
  )
  
  cat("\nSaved: ",path,"\n",sep="")
  
  p_final
}

# -----------------------------------------------------------------------
# 7. Fig. S6 and Fig. S7
# -----------------------------------------------------------------------
FigS6 <- make_hemi_fig(
  "Northern",
  "FigS6"
)

FigS7 <- make_hemi_fig(
  "Southern",
  "FigS7"
)

print(FigS6)
print(FigS7)