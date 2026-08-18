# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# Purpose: Fig. S4 - per-cell LDG slopes by geological era
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(purrr)
  library(palaeoverse); library(cowplot); library(patchwork); library(grid)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

set.seed(123)
occurrence_min <- 5
slope_method <- "balanced"

method_tag <- ifelse(
  slope_method=="balanced","per-cell balanced OLS","per-cell all-cells OLS"
)

slope_file <- sprintf(
  ifelse(
    slope_method=="balanced",
    "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS.csv",
    "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell all-cells OLS.csv"
  ),
  params$spacing,params$level,rich_params$n_lat_bins
)

dir.create("./figures",recursive=TRUE,showWarnings=FALSE)

hemi_cols <- c(
  Northern="#0072B2",
  Southern="#E69F00",
  `Poor quality`="#D3D3D3"
)
hemi_shapes <- c(Northern=16,Southern=17)

# -----------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------

rich_df <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing,params$level
))

LDG_slope <- read.csv(slope_file)
time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(n=rich_params$n_lat_bins) %>%
  arrange(min)

lat_zone_lookup <- lat_bins %>%
  mutate(
    lat_bin_mid=mid,
    abs_lat_bin_mid=round(abs(mid),6),
    lat_zone=case_when(
      abs_lat_bin_mid<30 ~ "tropical",
      abs_lat_bin_mid<60 ~ "temperate",
      abs_lat_bin_mid<=90 ~ "polar",
      TRUE ~ NA_character_
    )
  ) %>%
  select(lat_bin_mid,bin,abs_lat_bin_mid,lat_zone)

rich_df <- rich_df %>%
  filter(nT>=occurrence_min,t<=2*nT) %>%
  mutate(stage=time_bins$interval_name[match(bin_midpoint,time_bins$mid_ma)]) %>%
  filter(bin_midpoint<=486.85) %>%
  mutate(
    bin_index=findInterval(cell_lat,c(lat_bins$min,Inf)),
    bin=lat_bins$bin[bin_index],
    abs_lat=abs(cell_lat),
    hemisphere=ifelse(cell_lat>=0,"Northern","Southern")
  ) %>%
  left_join(lat_zone_lookup,by="bin") %>%
  filter(
    !is.na(bin),!is.na(abs_lat),!is.na(abs_lat_bin_mid),
    !is.na(lat_zone),!is.na(hemisphere),!is.na(stage)
  ) %>%
  group_by(bin_midpoint) %>%
  mutate(qD_normalized=qD*100/max(qD,na.rm=TRUE)) %>%
  ungroup()

# -----------------------------------------------------------------------
# 2. QC
# -----------------------------------------------------------------------

adjacent_df <- has_adjacent_bins(rich_df,lat_bins) %>%
  distinct(bin_midpoint,hemisphere,label) %>%
  transmute(
    bin_midpoint,hemisphere,
    has_adjacent_tt=label=="good"
  )

qc_label_df <- rich_df %>%
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
    hemisphere_mod=ifelse(label=="bad","Poor quality",hemisphere)
  )

rich_df <- rich_df %>%
  left_join(
    qc_label_df %>%
      select(bin_midpoint,hemisphere,label,hemisphere_mod),
    by=c("bin_midpoint","hemisphere")
  ) %>%
  mutate(
    label=coalesce(label,"bad"),
    hemisphere_mod=coalesce(hemisphere_mod,"Poor quality")
  )

# -----------------------------------------------------------------------
# 3. OLS lines
# -----------------------------------------------------------------------

LDG_slope <- LDG_slope %>%
  mutate(
    bin_midpoint=as.numeric(as.character(bin_midpoint)),
    slope=as.numeric(as.character(slope)),
    intercept=as.numeric(as.character(intercept)),
    hemisphere=as.character(hemisphere)
  ) %>%
  select(-any_of(c("label","color","hemisphere_mod"))) %>%
  left_join(
    qc_label_df %>%
      select(bin_midpoint,hemisphere,label,hemisphere_mod),
    by=c("bin_midpoint","hemisphere")
  ) %>%
  mutate(
    label=coalesce(label,"bad"),
    hemisphere_mod=coalesce(hemisphere_mod,"Poor quality")
  )

line_range <- rich_df %>%
  group_by(bin_midpoint,stage,hemisphere) %>%
  summarise(
    x_min=min(abs_lat,na.rm=TRUE),
    x_max=max(abs_lat,na.rm=TRUE),
    .groups="drop"
  )

ols_lines <- LDG_slope %>%
  left_join(line_range,by=c("bin_midpoint","stage","hemisphere")) %>%
  mutate(hemisphere_mod=ifelse(label=="bad","Poor quality",hemisphere)) %>%
  filter(!is.na(slope),!is.na(intercept),!is.na(x_min),!is.na(x_max)) %>%
  pmap_dfr(function(...){
    z <- tibble(...)
    x <- seq(z$x_min,z$x_max,length.out=100)
    tibble(
      bin_midpoint=z$bin_midpoint,
      stage=z$stage,
      hemisphere=z$hemisphere,
      hemisphere_mod=z$hemisphere_mod,
      abs_lat=x,
      fitted_values=z$intercept+z$slope*x
    )
  })

# -----------------------------------------------------------------------
# 4. Geological eras
# -----------------------------------------------------------------------

stage_order <- rich_df %>%
  distinct(bin_midpoint,stage) %>%
  arrange(desc(bin_midpoint)) %>%
  mutate(stage_id=row_number())

required <- c(
  "Tremadocian","Changhsingian",
  "Induan","Maastrichtian","Danian"
)

miss <- setdiff(required,stage_order$stage)
if(length(miss)>0)
  stop("Missing stage(s): ",paste(miss,collapse=", "))

i_trem  <- match("Tremadocian",stage_order$stage)
i_chang <- match("Changhsingian",stage_order$stage)
i_ind   <- match("Induan",stage_order$stage)
i_maas  <- match("Maastrichtian",stage_order$stage)
i_dan   <- match("Danian",stage_order$stage)

stage_order <- stage_order %>%
  mutate(
    era=case_when(
      between(stage_id,i_trem,i_chang) ~ "Palaeozoic",
      between(stage_id,i_ind,i_maas) ~ "Mesozoic",
      stage_id>=i_dan ~ "Cenozoic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(era))

rich_df <- rich_df %>%
  left_join(stage_order %>% select(bin_midpoint,era),by="bin_midpoint")

ols_lines <- ols_lines %>%
  left_join(stage_order %>% select(bin_midpoint,era),by="bin_midpoint")

# -----------------------------------------------------------------------
# 5. Era plot function
# -----------------------------------------------------------------------

make_era_plot <- function(era_name,show_x=FALSE){
  
  lev <- stage_order %>%
    filter(era==era_name) %>%
    arrange(desc(bin_midpoint)) %>%
    pull(stage)
  
  df <- rich_df %>%
    filter(era==era_name) %>%
    mutate(stage=factor(stage,levels=lev))
  
  ln <- ols_lines %>%
    filter(era==era_name) %>%
    mutate(stage=factor(stage,levels=lev))
  
  ggplot(
    df,
    aes(
      abs_lat,qD_normalized,
      color=hemisphere_mod,
      shape=hemisphere
    )
  )+
    geom_point(alpha=.65,size=1)+
    geom_line(
      data=ln,
      aes(
        abs_lat,fitted_values,
        color=hemisphere_mod,
        linetype=hemisphere_mod
      ),
      linewidth=.9,
      inherit.aes=FALSE
    )+
    scale_color_manual(
      values=hemi_cols,
      breaks=c("Northern","Southern","Poor quality")
    )+
    scale_shape_manual(
      values=hemi_shapes,
      breaks=c("Northern","Southern")
    )+
    scale_linetype_manual(
      values=c(
        Northern="solid",
        Southern="solid",
        `Poor quality`="solid"
      )
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
    facet_wrap(~stage,ncol=6,drop=FALSE)+
    labs(
      x=if(show_x) "Absolute palaeolatitude (°)" else NULL,
      y=NULL
    )+
    theme_minimal()+
    theme(
      strip.text=element_text(size=8,face="bold",margin=margin(1,1,1,1)),
      strip.placement="inside",
      panel.spacing.x=unit(.5,"lines"),
      panel.spacing.y=unit(.08,"lines"),
      panel.border=element_rect(color="black",fill=NA,linewidth=.8),
      panel.grid=element_blank(),
      axis.ticks=element_line(color="black",linewidth=.5),
      axis.ticks.length=unit(.06,"cm"),
      axis.text=element_text(size=8,color="black"),
      axis.title.x=element_text(size=12,color="black"),
      legend.position="none",
      
      # era dashed boundary
      plot.background=element_rect(
        color="grey40",
        fill=NA,
        linewidth=.9,
        linetype="dashed"
      ),
      
      # distance between facets and dashed boundary
      plot.margin=margin(7,7,7,7)
    )
}

era_label <- function(x){
  ggplot()+
    annotate(
      "text",x=.5,y=.5,label=x,
      angle=270,fontface="bold",size=3.5
    )+
    xlim(0,1)+ylim(0,1)+theme_void()
}

# -----------------------------------------------------------------------
# 6. Three era blocks
# -----------------------------------------------------------------------

p_pal <- make_era_plot("Palaeozoic",FALSE)
p_mes <- make_era_plot("Mesozoic",FALSE)
p_cen <- make_era_plot("Cenozoic",TRUE)

row_pal <- p_pal + era_label("Palaeozoic") +
  plot_layout(widths=c(1,.045))

row_mes <- p_mes + era_label("Mesozoic") +
  plot_layout(widths=c(1,.045))

row_cen <- p_cen + era_label("Cenozoic") +
  plot_layout(widths=c(1,.045))

# -----------------------------------------------------------------------
# 7. Legend
# -----------------------------------------------------------------------

legend_plot <- ggplot(
  rich_df,
  aes(
    abs_lat,qD_normalized,
    color=hemisphere_mod,
    shape=hemisphere
  )
)+
  geom_point()+
  geom_line(aes(group=hemisphere_mod,linetype=hemisphere_mod))+
  scale_color_manual(
    name="Line",
    values=hemi_cols,
    breaks=c("Northern","Southern","Poor quality")
  )+
  scale_shape_manual(
    name="Cell",
    values=hemi_shapes,
    breaks=c("Northern","Southern")
  )+
  scale_linetype_manual(
    values=c(
      Northern="solid",
      Southern="solid",
      `Poor quality`="solid"
    )
  )+
  guides(
    shape=guide_legend(
      title="Cell",
      override.aes=list(color="black",size=2)
    ),
    color=guide_legend(
      title="Line",
      override.aes=list(shape=NA,linewidth=1)
    ),
    linetype="none"
  )+
  theme_void()+
  theme(
    legend.position="bottom",
    legend.box="horizontal",
    legend.title=element_text(size=9),
    legend.text=element_text(size=8)
  )

legend_patch <- wrap_elements(
  full=cowplot::get_legend(legend_plot)
)

# -----------------------------------------------------------------------
# 8. Combine Fig. S4
# -----------------------------------------------------------------------

gap <- plot_spacer()

n_pal <- ceiling(sum(stage_order$era=="Palaeozoic")/6)
n_mes <- ceiling(sum(stage_order$era=="Mesozoic")/6)
n_cen <- ceiling(sum(stage_order$era=="Cenozoic")/6)

FigS4 <-
  row_pal /
  gap /
  row_mes /
  gap /
  row_cen /
  legend_patch +
  plot_layout(
    heights=c(
      n_pal,
      .12,
      n_mes,
      .12,
      n_cen,
      .45
    )
  )

FigS4 <- (
  wrap_elements(
    full=textGrob(
      "Normalized generic richness",
      rot=90,
      gp=gpar(fontsize=12)
    )
  ) |
    FigS4
)+
  plot_layout(widths=c(.035,1))

# -----------------------------------------------------------------------
# 9. Save ONLY Fig. S4
# -----------------------------------------------------------------------

FigS4_path <- sprintf(
  "./figures/FigS4_%skm_%squota_%s_equal-area_latbins_%s.jpg",
  params$spacing,
  params$level,
  rich_params$n_lat_bins,
  slope_method
)

if(file.exists(FigS4_path))
  file.remove(FigS4_path)

ggsave(
  filename=FigS4_path,
  plot=FigS4,
  width=8,
  height=9,
  dpi=300,
  bg="white"
)

print(FigS4)
cat("\nSaved:",normalizePath(FigS4_path),"\n")