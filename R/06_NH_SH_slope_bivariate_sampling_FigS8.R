# -----------------------------------------------------------------------
# Project: LDG_climate_state
# File: 06_NH_SH_slope_bivariate_sampling_FigS8.R
# Purpose: NH-SH signed slope comparison and sampling diagnostics
#          for per-cell balanced resampling LDG slopes
# -----------------------------------------------------------------------

# rm(list=ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(deeptime)
  library(grid)
})

source("./R/options.R")

# Disable showtext if it was enabled elsewhere
if ("showtext" %in% loadedNamespaces()) {
  showtext::showtext_auto(FALSE)
}

era_uk <- deeptime::get_scale_data("era")
era_uk$name[era_uk$name=="Paleozoic"] <- "Palaeozoic"

# -----------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------

baseline_qc <- "occurrence5_k1_tropical_temperate"
method_use <- "per_cell_balanced_resampling_OLS"
metric_use <- "median_resampled_slope"

era_cols <- c(
  "Palaeozoic"="#9BBA7F",
  "Mesozoic"="#67C5CA",
  "Cenozoic"="#F2D2A2"
)

percentile_tag <- ifelse(
  length(rich_params$percentiles)>1,
  "allq",
  rich_params$percentiles[1]
)

analysis_tag <- sprintf(
  "%skm_%squota_%slat_%s",
  params$spacing,
  params$level,
  rich_params$n_lat_bins,
  percentile_tag
)

out_dir <- "./results/NH_SH_slope_diagnostics_percell_balanced"
fig_jpg_dir <- "./figures/jpg/NH_SH_slope_diagnostics_percell_balanced"
fig_pdf_dir <- "./figures/pdf/NH_SH_slope_diagnostics_percell_balanced"

sampling_ts_jpg_dir <- file.path(
  fig_jpg_dir,
  "sampling_profile_dissim_time_series"
)

sampling_ts_pdf_dir <- file.path(
  fig_pdf_dir,
  "sampling_profile_dissim_time_series"
)

dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(fig_jpg_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(fig_pdf_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(sampling_ts_jpg_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(sampling_ts_pdf_dir,recursive=TRUE,showWarnings=FALSE)

# -----------------------------------------------------------------------
# 2. Read data
# -----------------------------------------------------------------------

slope_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS.csv",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

if(!file.exists(slope_path)){
  stop(
    "Cannot find baseline slope file:\n",
    slope_path
  )
}

LDG_slope <- read.csv(slope_path)
time_bins <- readRDS("./data/time_bins.RDS")

meta_df <- time_bins %>%
  mutate(
    bin_midpoint=mid_ma,
    stage=interval_name,
    era=case_when(
      bin_midpoint<66~"Cenozoic",
      bin_midpoint>=66 & bin_midpoint<251.902~"Mesozoic",
      bin_midpoint>=251.902~"Palaeozoic",
      TRUE~NA_character_
    )
  ) %>%
  select(bin_midpoint,stage,era)

# -----------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------

cor_fun <- function(x,y){
  
  ok <- complete.cases(x,y)
  
  if(sum(ok)<3){
    return(tibble(
      n=sum(ok),
      estimate=NA_real_,
      p_value=NA_real_
    ))
  }
  
  z <- suppressWarnings(
    cor.test(
      x[ok],y[ok],
      method="spearman",
      exact=FALSE
    )
  )
  
  tibble(
    n=sum(ok),
    estimate=unname(z$estimate),
    p_value=z$p.value
  )
}

cor_one <- function(data,x,y,test_name){
  
  cor_fun(data[[x]],data[[y]]) %>%
    mutate(
      test=test_name,
      x_var=x,
      y_var=y
    ) %>%
    select(
      test,x_var,y_var,
      n,estimate,p_value
    )
}

fmt_num <- function(x,digits=3){
  ifelse(
    is.na(x),
    "NA",
    format(
      round(x,digits),
      nsmall=digits
    )
  )
}

safe_range_expand <- function(x,expand_frac=.08){
  
  x <- x[is.finite(x)]
  
  if(!length(x))
    return(c(-1,1))
  
  xr <- range(x)
  
  if(xr[1]==xr[2]){
    xr+c(-.5,.5)
  }else{
    xr+c(-1,1)*diff(xr)*expand_frac
  }
}

# -----------------------------------------------------------------------
# Themes: Arial, without showtext
# -----------------------------------------------------------------------

base_theme <- theme_minimal(base_family="Arial") +
  theme(
    text=element_text(family="Arial"),
    panel.grid=element_blank(),
    panel.border=element_rect(
      colour="black",
      fill=NA,
      linewidth=.6
    ),
    axis.text=element_text(size=10,colour="black"),
    axis.title=element_text(size=12,colour="black"),
    axis.ticks=element_line(colour="black"),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.title=element_text(size=11,face="bold"),
    legend.text=element_text(size=10),
    plot.tag=element_text(size=14,face="bold"),
    aspect.ratio=1,
    
    # reduce outer margins
    plot.margin=margin(2,2,2,2)
  )

ts_theme <- theme_minimal(
  base_family="Arial"
)+
  theme(
    text=element_text(family="Arial"),
    panel.grid=element_blank(),
    panel.border=element_rect(
      color="black",
      fill=NA,
      linewidth=.6
    ),
    axis.title=element_text(
      size=12,
      colour="black"
    ),
    axis.text=element_text(
      size=10,
      colour="black"
    ),
    axis.ticks=element_line(
      color="black",
      linewidth=.5
    ),
    legend.position="bottom",
    legend.title=element_text(
      size=11,
      face="bold"
    ),
    legend.text=element_text(
      size=10
    ),
    plot.title=element_text(
      size=12,
      face="bold",
      hjust=.5
    )
  )

# -----------------------------------------------------------------------
# 4. Sampling-profile dissimilarity
# -----------------------------------------------------------------------

sampling_dissim <- LDG_slope %>%
  mutate(
    bin_midpoint=as.numeric(bin_midpoint),
    sampling_profile_dissim=as.numeric(
      sampling_profile_dissim
    )
  ) %>%
  filter(
    method_group==method_use,
    slope_metric==metric_use,
    qc_name==baseline_qc
  ) %>%
  select(
    bin_midpoint,
    sampling_profile_dissim
  ) %>%
  distinct(
    bin_midpoint,
    .keep_all=TRUE
  )

write.csv(
  sampling_dissim,
  file.path(
    out_dir,
    paste0(
      analysis_tag,
      "_sampling_profile_dissim_",
      baseline_qc,
      ".csv"
    )
  ),
  row.names=FALSE
)

# -----------------------------------------------------------------------
# 5. Paired NH-SH data
# -----------------------------------------------------------------------

paired_df <- LDG_slope %>%
  mutate(
    bin_midpoint=as.numeric(bin_midpoint),
    slope=as.numeric(slope),
    k_median=as.numeric(k_median),
    k_cv=as.numeric(k_cv),
    k_mean=as.numeric(k_mean)
  ) %>%
  filter(
    method_group==method_use,
    slope_metric==metric_use,
    qc_name==baseline_qc,
    label=="good",
    hemisphere %in% c(
      "Northern",
      "Southern"
    )
  ) %>%
  left_join(
    meta_df,
    by=c("bin_midpoint","stage")
  ) %>%
  select(
    bin_midpoint,
    stage,
    era,
    hemisphere,
    slope,
    k_median,
    k_cv,
    k_mean
  ) %>%
  pivot_wider(
    names_from=hemisphere,
    values_from=c(
      slope,
      k_median,
      k_cv,
      k_mean
    ),
    names_sep="_"
  ) %>%
  distinct(
    bin_midpoint,
    .keep_all=TRUE
  ) %>%
  left_join(
    sampling_dissim,
    by="bin_midpoint"
  ) %>%
  filter(
    !is.na(slope_Northern),
    !is.na(slope_Southern)
  ) %>%
  mutate(
    slope_diff_abs=
      abs(slope_Northern-slope_Southern),
    
    slope_diff_signed=
      slope_Northern-slope_Southern,
    
    k_cv_diff=
      abs(k_cv_Northern-k_cv_Southern),
    
    k_cv_diff_signed=
      k_cv_Northern-k_cv_Southern,
    
    k_cv_asym=
      abs(log(
        (k_cv_Northern+1)/
          (k_cv_Southern+1)
      )),
    
    k_median_diff=
      abs(
        k_median_Northern-
          k_median_Southern
      ),
    
    k_median_asym=
      abs(log(
        (k_median_Northern+1)/
          (k_median_Southern+1)
      )),
    
    k_median_diff_signed=
      k_median_Northern-
      k_median_Southern,
    
    k_mean_diff=
      abs(
        k_mean_Northern-
          k_mean_Southern
      ),
    
    k_mean_asym=
      abs(log(
        (k_mean_Northern+1)/
          (k_mean_Southern+1)
      )),
    
    k_mean_diff_signed=
      k_mean_Northern-
      k_mean_Southern,
    
    era=factor(
      era,
      levels=c(
        "Palaeozoic",
        "Mesozoic",
        "Cenozoic"
      )
    ),
    
    method_group=method_use,
    slope_metric=metric_use,
    qc_name=baseline_qc
  )

write.csv(
  paired_df,
  file.path(
    out_dir,
    paste0(
      analysis_tag,
      "_paired_NH_SH_percell_balanced_",
      baseline_qc,
      ".csv"
    )
  ),
  row.names=FALSE
)

if(nrow(paired_df)<3){
  stop(
    "Fewer than 3 paired NH-SH rows. ",
    "Cannot run diagnostics."
  )
}

# -----------------------------------------------------------------------
# 6. Correlations
# -----------------------------------------------------------------------

cor_results <- bind_rows(
  
  cor_one(
    paired_df,
    "slope_Northern",
    "slope_Southern",
    "NH slope vs SH slope"
  ),
  
  cor_one(
    paired_df,
    "sampling_profile_dissim",
    "slope_diff_abs",
    "|NH-SH slope| vs sampling-profile dissimilarity"
  ),
  
  cor_one(
    paired_df,
    "k_cv_diff",
    "slope_diff_abs",
    "|NH-SH slope| vs k_cv difference"
  ),
  
  cor_one(
    paired_df,
    "k_median_diff",
    "slope_diff_abs",
    "|NH-SH slope| vs absolute k_median difference"
  ),
  
  cor_one(
    paired_df,
    "k_median_asym",
    "slope_diff_abs",
    "|NH-SH slope| vs k_median asymmetry"
  ),
  
  cor_one(
    paired_df,
    "sampling_profile_dissim",
    "slope_diff_signed",
    "NH-SH slope vs sampling-profile dissimilarity"
  ),
  
  cor_one(
    paired_df,
    "k_cv_diff_signed",
    "slope_diff_signed",
    "NH-SH slope vs signed k_cv difference"
  ),
  
  cor_one(
    paired_df,
    "k_median_diff_signed",
    "slope_diff_signed",
    "NH-SH slope vs signed k_median difference"
  )
  
) %>%
  mutate(
    method_group=method_use,
    slope_metric=metric_use,
    qc_name=baseline_qc,
    estimate_label=fmt_num(estimate),
    p_value_label=fmt_num(p_value)
  ) %>%
  select(
    method_group,
    slope_metric,
    qc_name,
    everything()
  )

write.csv(
  cor_results,
  file.path(
    out_dir,
    paste0(
      analysis_tag,
      "_NH_SH_correlations_percell_balanced_",
      baseline_qc,
      ".csv"
    )
  ),
  row.names=FALSE
)

print(cor_results)

# -----------------------------------------------------------------------
# 7. Sampling-profile dissimilarity time series
# -----------------------------------------------------------------------

draw_sampling_dissim_time_series <- function(sampling_df){
  
  message(
    "Drawing sampling-profile dissimilarity time series: ",
    baseline_qc
  )
  
  time_bins_use <- time_bins %>%
    filter(
      mid_ma<=486.85,
      mid_ma>=0
    )
  
  sampling_ts_df <- sampling_df %>%
    mutate(
      bin_midpoint=
        as.numeric(bin_midpoint)
    ) %>%
    select(
      bin_midpoint,
      sampling_profile_dissim
    ) %>%
    distinct() %>%
    right_join(
      time_bins_use %>%
        transmute(
          bin_midpoint=mid_ma,
          stage=interval_name,
          era=case_when(
            bin_midpoint<66~
              "Cenozoic",
            bin_midpoint>=66 &
              bin_midpoint<251.902~
              "Mesozoic",
            TRUE~
              "Palaeozoic"
          )
        ),
      by="bin_midpoint"
    ) %>%
    mutate(
      era=factor(
        era,
        levels=c(
          "Palaeozoic",
          "Mesozoic",
          "Cenozoic"
        )
      )
    ) %>%
    arrange(desc(bin_midpoint))
  
  if(
    sum(
      !is.na(
        sampling_ts_df$
        sampling_profile_dissim
      )
    )<3
  ){
    message(
      "Skip sampling-profile dissimilarity time series: ",
      "fewer than 3 non-NA rows."
    )
    return(NULL)
  }
  
  data(periods)
  data(epochs)
  
  x_max_val <- max(
    time_bins_use$max_ma,
    na.rm=TRUE
  )
  
  y_max_val <- max(
    sampling_ts_df$
      sampling_profile_dissim,
    na.rm=TRUE
  )
  
  if(
    !is.finite(y_max_val) ||
    y_max_val<=0
  ){
    y_max_val <- 1
  }
  
  major_boundaries <- periods$max_age
  
  major_boundaries <- major_boundaries[
    is.finite(major_boundaries) &
      major_boundaries>=0 &
      major_boundaries<=x_max_val
  ]
  
  p <- ggplot(
    sampling_ts_df,
    aes(
      x=bin_midpoint,
      y=sampling_profile_dissim
    )
  )+
    geom_vline(
      xintercept=major_boundaries,
      color="black",
      linewidth=.35,
      alpha=.75
    )+
    geom_line(
      linewidth=.8,
      color="black",
      na.rm=FALSE
    )+
    geom_point(
      aes(fill=era),
      shape=21,
      size=2.4,
      stroke=.45,
      color="black",
      na.rm=TRUE
    )+
    scale_fill_manual(
      values=era_cols,
      drop=FALSE,
      name="Era"
    )+
    scale_x_reverse(
      limits=c(x_max_val,0),
      breaks=seq(450,0,-50),
      expand=c(0,0)
    )+
    scale_y_continuous(
      limits=c(0,y_max_val),
      breaks=pretty(
        c(0,y_max_val),
        n=5
      ),
      expand=c(0,0)
    )+
    coord_geo(
      xlim=c(x_max_val,0),
      pos=as.list(rep("bottom",2)),
      dat=list("periods",era_uk),
      height=list(
        unit(1.35,"lines"),
        unit(1.35,"lines")
      ),
      lab_color="black",
      rot=list(0,0),
      abbrv=list(TRUE,FALSE)
    )+
    labs(
      x="Time (Ma)",
      y="Sampling-profile dissimilarity"
    )+
    ts_theme+
    theme(
      legend.box.margin=margin(
        t=-6,r=0,b=0,l=0
      ),
      legend.margin=margin(
        t=-8,r=0,b=0,l=0
      ),
      legend.spacing.y=unit(
        .05,"cm"
      ),
      plot.margin=margin(
        5,5,2,5
      )
    )
  
  print(p)
  
  jpg_path <- file.path(
    sampling_ts_jpg_dir,
    paste0(
      analysis_tag,
      "_sampling_profile_dissim_time_series_",
      "percell_balanced_",
      baseline_qc,
      ".jpg"
    )
  )
  
  pdf_path <- file.path(
    sampling_ts_pdf_dir,
    paste0(
      analysis_tag,
      "_sampling_profile_dissim_time_series_",
      "percell_balanced_",
      baseline_qc,
      ".pdf"
    )
  )
  
  ggsave(
    filename=jpg_path,
    plot=p,
    width=8,
    height=4,
    dpi=300
  )
  
  ggsave(
    filename=pdf_path,
    plot=p,
    width=8,
    height=4,
    dpi=300,
    device=cairo_pdf
  )
  
  print(jpg_path)
  print(pdf_path)
  
  p
}

sampling_ts_plot <-
  draw_sampling_dissim_time_series(
    sampling_dissim
  )

# -----------------------------------------------------------------------
# 8. Bivariate diagnostic plot
# -----------------------------------------------------------------------

cor_A <- cor_results %>% filter(test=="NH slope vs SH slope")
cor_B <- cor_results %>% filter(
  test=="|NH-SH slope| vs sampling-profile dissimilarity"
)

axis_lim <- safe_range_expand(
  c(paired_df$slope_Northern,paired_df$slope_Southern)
)

# -----------------------------------------------------------------------
# p1 annotation: bottom-right, right aligned
# -----------------------------------------------------------------------
xA <- axis_lim[2]-.06*diff(axis_lim)
yA <- axis_lim[1]+.07*diff(axis_lim)
dyA <- .075*diff(axis_lim)

lab_A_df <- tibble(
  x=xA,
  y=c(yA+2*dyA,yA+dyA,yA),
  label=c(
    paste0("\u03c1 = ",cor_A$estimate_label),
    paste0("p = ",cor_A$p_value_label),
    paste0("n = ",cor_A$n)
  )
)

p1 <- ggplot(
  paired_df,
  aes(slope_Northern,slope_Southern,fill=era)
)+
  geom_hline(
    yintercept=0,
    linetype="dashed",
    colour="grey40",
    linewidth=.4
  )+
  geom_vline(
    xintercept=0,
    linetype="dashed",
    colour="grey40",
    linewidth=.4
  )+
  geom_abline(
    slope=1,
    intercept=0,
    linetype="dashed",
    linewidth=.6
  )+
  geom_point(
    shape=21,size=3,stroke=.6,
    colour="black",alpha=.9
  )+
  geom_text(
    data=lab_A_df,
    aes(x=x,y=y,label=label),
    inherit.aes=FALSE,
    hjust=1,
    size=3.5,
    colour="black"
  )+
  scale_fill_manual(
    values=era_cols,
    drop=FALSE
  )+
  coord_equal(
    xlim=axis_lim,
    ylim=axis_lim,
    expand=FALSE
  )+
  labs(
    x="Northern Hemisphere slope",
    y="Southern Hemisphere slope",
    fill="Era",
    tag="(a)"
  )+
  base_theme

# -----------------------------------------------------------------------
# p2
# -----------------------------------------------------------------------
fit_ok <- complete.cases(
  paired_df$sampling_profile_dissim,
  paired_df$slope_diff_abs
)

ok_fit <- sum(fit_ok)>=3 &&
  length(unique(paired_df$sampling_profile_dissim[fit_ok]))>1 &&
  length(unique(paired_df$slope_diff_abs[fit_ok]))>1

p2 <- ggplot(
  paired_df,
  aes(sampling_profile_dissim,slope_diff_abs,fill=era)
)+
  geom_point(
    shape=21,size=3,stroke=.6,
    colour="black",alpha=.9
  )+
  scale_fill_manual(
    values=era_cols,
    drop=FALSE
  )+
  labs(
    x="NH-SH sampling-profile dissimilarity",
    y="|Northern slope - Southern slope|",
    fill="Era",
    tag="(b)"
  )+
  base_theme

if(ok_fit){
  p2 <- p2+
    geom_smooth(
      aes(group=1),
      method="lm",
      formula=y~x,
      se=TRUE,
      colour="black",
      fill="grey80",
      linewidth=.6
    )
}

# -----------------------------------------------------------------------
# p2 annotation: upper-left region, right aligned
# -----------------------------------------------------------------------
xrB <- range(
  paired_df$sampling_profile_dissim,
  na.rm=TRUE
)

yrB <- range(
  paired_df$slope_diff_abs,
  na.rm=TRUE
)

xB <- xrB[1]+.26*diff(xrB)
yB <- yrB[2]-.03*diff(yrB)
dyB <- .075*diff(yrB)

lab_B_df <- tibble(
  x=xB,
  y=c(yB,yB-dyB,yB-2*dyB),
  label=c(
    paste0("\u03c1 = ",cor_B$estimate_label),
    paste0("p = ",cor_B$p_value_label),
    paste0("n = ",cor_B$n)
  )
)

p2 <- p2+
  geom_text(
    data=lab_B_df,
    aes(x=x,y=y,label=label),
    inherit.aes=FALSE,
    hjust=1,
    size=3.5,
    colour="black"
  )

# -----------------------------------------------------------------------
# Final plot
# -----------------------------------------------------------------------
final_plot <- wrap_plots(
  p1,p2,
  ncol=2,
  guides="collect"
)&
  theme(
    legend.position="bottom",
    legend.direction="horizontal"
  )

print(final_plot)

jpg_final <- file.path(
  fig_jpg_dir,
  paste0(
    analysis_tag,
    "_NH_SH_sampling_diagnostics_percell_balanced_",
    baseline_qc,
    ".jpg"
  )
)

pdf_final <- file.path(
  fig_pdf_dir,
  paste0(
    analysis_tag,
    "_NH_SH_sampling_diagnostics_percell_balanced_",
    baseline_qc,
    ".pdf"
  )
)

ggsave(
  jpg_final,
  final_plot,
  width=7,
  height=4,
  dpi=300,
  bg="white"
)

ggsave(
  pdf_final,
  final_plot,
  width=7,
  height=4,
  dpi=300,
  device=cairo_pdf,
  bg="white"
)

print(jpg_final)
print(pdf_final)
# -----------------------------------------------------------------------
# 9. Run summary
# -----------------------------------------------------------------------

run_summary <- tibble(
  method_group=method_use,
  slope_metric=metric_use,
  qc_name=baseline_qc,
  n_pairs=nrow(paired_df),
  status="ok"
)

write.csv(
  run_summary,
  file.path(
    out_dir,
    paste0(
      analysis_tag,
      "_NH_SH_run_summary_percell_balanced_",
      baseline_qc,
      ".csv"
    )
  ),
  row.names=FALSE
)

message(
  "NH-SH per-cell balanced slope diagnostics finished."
)