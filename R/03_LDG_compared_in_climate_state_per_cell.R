# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 03_LDG_compared_in_climate_state_per_cell_balanced.R
# Purpose: Read per-cell balanced resampling LDG slopes, compare among
#          climate states, and draw time-series + boxplot
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(grid)
  library(patchwork)
  library(cowplot)
  library(deeptime)
})
era_uk <- deeptime::get_scale_data("era")
era_uk$name[era_uk$name == "Paleozoic"] <- "Palaeozoic"
source("./R/options.R")

# -----------------------------------------------------------------------
# 0.Settings
# -----------------------------------------------------------------------

climate_levels<-c("Coldhouse","Coolhouse","Transitional","Warmhouse","Hothouse")
# climate_colors<-c(
#   "Coldhouse"="#005344",
#   "Coolhouse"="#007d65",
#   "Transitional"="#c8c7c7",
#   "Warmhouse"="#b57a51",
#   "Hothouse"="#95484b"
# )
climate_colors <- c(
  "Coldhouse"    = "#B2CCC7",
  "Coolhouse"    = "#B3D8D0",
  "Transitional" = "#E4E4E4",
  "Warmhouse"    = "#DECABC",
  "Hothouse"     = "#D6B6B7"
)
hemi_cols<-c("Northern"="#0072B2","Southern"="#E69F00")

dir.create("./results",recursive=TRUE,showWarnings=FALSE)
dir.create("./figures/jpg",recursive=TRUE,showWarnings=FALSE)
dir.create("./figures/pdf",recursive=TRUE,showWarnings=FALSE)

# -----------------------------------------------------------------------
# 1.Helper functions
# -----------------------------------------------------------------------

safe_shapiro <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(unique(x)) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}

safe_wilcox_less <- function(x, y, min_group_n = 2) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  if (length(x) < min_group_n || length(y) < min_group_n) {
    return(list(
      p_value = NA_real_,
      statistic = NA_real_,
      test_note = paste0("Skipped: n < ", min_group_n)
    ))
  }
  
  wt <- suppressWarnings(
    wilcox.test(x, y, alternative = "less", exact = FALSE)
  )
  
  list(
    p_value = wt$p.value,
    statistic = as.numeric(wt$statistic),
    test_note = "Tested"
  )
}

safe_num <- function(x) {
  x <- as.numeric(x)
  x[is.nan(x) | is.infinite(x)] <- NA_real_
  x
}

# -----------------------------------------------------------------------
# 2.Read data
# -----------------------------------------------------------------------

LDG_slope<-read.csv(sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS.csv",
  params$spacing,params$level,rich_params$n_lat_bins
))

climate_states<-read.csv("./data/climate_states.csv")
time_bins<-readRDS("./data/time_bins.RDS")

climate_states<-climate_states%>%
  filter(bottom<=486.8500,top>=0)%>%
  mutate(
    climate_state=factor(climate_state,levels=climate_levels),
    climate_color=unname(climate_colors[as.character(climate_state)])
  )

time_bins<-time_bins%>%filter(min_ma<=486.8500,max_ma>=0)

# -----------------------------------------------------------------------
# 3.Merge with climate states
# -----------------------------------------------------------------------

slope_cli_df <- LDG_slope %>%
  mutate(
    bin_midpoint = safe_num(bin_midpoint),
    slope = safe_num(slope),
    slope_lower_95 = safe_num(slope_lower_95),
    slope_upper_95 = safe_num(slope_upper_95)
  ) %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid"))

# Compatibility with older baseline output
if (!"qc_name" %in% names(slope_cli_df)) {
  slope_cli_df$qc_name <- "occurrence5_k1_tropical_temperate"
}

if (!"method_group" %in% names(slope_cli_df)) {
  slope_cli_df$method_group <- "per_cell_balanced_resampling_OLS"
}

if (!"slope_metric" %in% names(slope_cli_df)) {
  slope_cli_df$slope_metric <- "median_resampled_slope"
}

slope_cli_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins LDG slope per-cell balanced OLS and climate states.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(slope_cli_df, slope_cli_path, row.names = FALSE)

slope_cli_df_filter <- slope_cli_df %>%
  filter(
    qc_name == "occurrence5_k1_tropical_temperate",
    method_group == "per_cell_balanced_resampling_OLS",
    slope_metric == "median_resampled_slope",
    climate_state %in% climate_levels,
    label == "good",
    !is.na(slope)
  )
# -----------------------------------------------------------------------
# 4.Climate-state tests
# -----------------------------------------------------------------------

climate_pairs <- combn(climate_levels, 2, simplify = FALSE)

test_results <- map_df(climate_pairs, function(pair) {
  
  g1 <- slope_cli_df_filter %>%
    filter(climate_state == pair[1]) %>%
    pull(slope) %>%
    as.numeric()
  
  g2 <- slope_cli_df_filter %>%
    filter(climate_state == pair[2]) %>%
    pull(slope) %>%
    as.numeric()
  
  shapiro1_p <- safe_shapiro(g1)
  shapiro2_p <- safe_shapiro(g2)
  
  normal <- !is.na(shapiro1_p) &&
    !is.na(shapiro2_p) &&
    shapiro1_p > 0.05 &&
    shapiro2_p > 0.05
  
  tibble(
    method_group = "per_cell_balanced_resampling_OLS",
    slope_metric = "median_resampled_slope",
    qc_name = "occurrence5_k1_tropical_temperate",
    group1 = pair[1],
    group2 = pair[2],
    normal = normal,
    test_type = ifelse(normal, "t-test", "Wilcoxon"),
    shapiro_p1 = shapiro1_p,
    shapiro_p2 = shapiro2_p
  )
})

test_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins per-cell balanced OLS normality tests.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(test_results, test_path, row.names = FALSE)

wil_results <- map_df(climate_pairs, function(pair) {
  
  g1 <- slope_cli_df_filter %>%
    filter(climate_state == pair[1]) %>%
    pull(slope) %>%
    as.numeric()
  
  g2 <- slope_cli_df_filter %>%
    filter(climate_state == pair[2]) %>%
    pull(slope) %>%
    as.numeric()
  
  g1 <- g1[!is.na(g1)]
  g2 <- g2[!is.na(g2)]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  wt <- safe_wilcox_less(g1, g2, min_group_n = 2)
  
  med1 <- ifelse(n1 > 0, median(g1, na.rm = TRUE), NA_real_)
  med2 <- ifelse(n2 > 0, median(g2, na.rm = TRUE), NA_real_)
  
  tibble(
    method_group = "per_cell_balanced_resampling_OLS",
    slope_metric = "median_resampled_slope",
    qc_name = "occurrence5_k1_tropical_temperate",
    group1 = pair[1],
    group2 = pair[2],
    n1 = n1,
    n2 = n2,
    median1 = med1,
    median2 = med2,
    median_diff = med1 - med2,
    p_value = wt$p_value,
    w_statistic = wt$statistic,
    iqr1 = ifelse(n1 > 0, IQR(g1, na.rm = TRUE), NA_real_),
    iqr2 = ifelse(n2 > 0, IQR(g2, na.rm = TRUE), NA_real_),
    test_note = wt$test_note
  )
}) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adjusted)

wil_path <- sprintf(
  "./results/%skm %squota %s equal-area latitude bins per-cell balanced OLS wilcoxon tests.csv",
  params$spacing, params$level, rich_params$n_lat_bins
)

write.csv(wil_results, wil_path, row.names = FALSE)

print(wil_results)
# -----------------------------------------------------------------------
# 5.Time-series data
# -----------------------------------------------------------------------
slope_data<-slope_cli_df%>%
  mutate(
    slope_type=case_when(
      hemisphere=="Northern"~"Northern",
      hemisphere=="Southern"~"Southern",
      TRUE~NA_character_
    ),
    uncertain = case_when(
      label == "bad" ~ NA,
      is.na(slope_lower_95) | is.na(slope_upper_95) ~ NA,
      slope_lower_95 <= 0 & slope_upper_95 >= 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    slope_value=ifelse(label=="bad",NA_real_,slope),
    slope_lower_plot=ifelse(label=="bad",NA_real_,slope_lower_95),
    slope_upper_plot=ifelse(label=="bad",NA_real_,slope_upper_95)
  )

data(periods)
data(epochs)

major_boundaries<-periods$max_age
x_max_val<-max(time_bins$max_ma,na.rm=TRUE)
x_min_val<-min(time_bins$min_ma,na.rm=TRUE)
major_boundaries_plot<-major_boundaries[major_boundaries>=0&major_boundaries<=x_max_val]
# 
y_min_val<-min(c(slope_data$slope_value,slope_data$slope_lower_plot),na.rm=TRUE)*1.01
y_max_val<-max(c(slope_data$slope_value,slope_data$slope_upper_plot),na.rm=TRUE)*1.1

# y_min_val<-min(c(slope_data$slope_value),na.rm=TRUE)*1
# y_max_val<-max(c(slope_data$slope_value),na.rm=TRUE)*1.1
if(!is.finite(y_min_val)||!is.finite(y_max_val)||y_min_val==y_max_val){
  y_min_val<--1
  y_max_val<-1
}

y_breaks<-pretty(c(y_min_val,y_max_val),n=5)

climate_shade_layer<-list(
  geom_rect(
    data=climate_states,
    aes(xmin=bottom,xmax=top,ymin=y_min_val,ymax=y_max_val*1.1,fill=I(climate_color)),
    inherit.aes=FALSE,colour=NA #,alpha=0.55
  )
)

plot_slope_ts<-function(df,hemi,col,tag_lab,show_top_bar=FALSE,show_x_tick=FALSE,bottom_space=0){
  p<-ggplot(filter(df,slope_type==hemi),aes(x=bin_midpoint,y=slope_value))+
    climate_shade_layer+
    geom_vline(xintercept=major_boundaries_plot,color="black",linewidth=0.35,alpha=0.7)
  
  if(show_top_bar){
    p<-p+
      geom_rect(
        data=climate_states,
        aes(xmin=bottom,xmax=top,ymin=y_max_val,ymax=y_max_val*1.1),
        fill=I(climate_states$climate_color),
        color="black",linewidth=0.3,inherit.aes=FALSE
      )
  }
  
  p<-p+
    geom_hline(yintercept=0,color="black",linewidth=0.4,linetype="dashed")+
    # geom_ribbon(
    #   aes(ymin=slope_lower_plot,ymax=slope_upper_plot),
    #   fill=unname(col),alpha=0.55,colour=NA,na.rm=TRUE
    # )+
    geom_errorbar(
      aes(ymin = slope_lower_plot, ymax = slope_upper_plot),
      width = 2, linewidth = 0.75, color = unname(col), na.rm = TRUE) +
    # geom_linerange(
    #   aes(ymin = slope_lower_plot, ymax = slope_upper_plot),
    #   colour = unname(col), alpha = 0.8, width = 2, linewidth = 0.75, na.rm = TRUE
    # ) +
    geom_line(linewidth=1,color=unname(col),na.rm=TRUE)+
    geom_point(
      aes(shape = abs(slope_value) < 0.1), #aes(shape=uncertain),
      size=2,stroke=0.5,color="black",fill=unname(col),na.rm=TRUE
    )+
    scale_shape_manual(values=c(`TRUE`=2,`FALSE`=21),na.translate=FALSE)+
    annotate(
      "rect",
      xmin=x_min_val,xmax=x_max_val, ymin=y_min_val,ymax=y_max_val*1.1,
      fill=NA,color="black",linewidth=0.8
    )+
    scale_x_reverse(
      name=NULL,
      limits=c(x_max_val,0),
      breaks=seq(500,0,-50),
      expand=c(0,0)
    )+
    scale_y_continuous(
      limits=c(y_min_val,y_max_val*1.1),
      breaks=y_breaks,
      expand=c(0,0)
    )+
    labs(
      x=NULL,
      y="Slope value",
      tag=tag_lab,
      subtitle=paste(hemi,"Hemisphere")
    )+
    theme_minimal()+
    theme(
      panel.grid=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=14),
      axis.text.y=element_text(size=14),
      axis.ticks.y=element_line(color="black",linewidth=0.5),
      plot.subtitle=element_text(size=12,face="bold"),
      plot.tag=element_text(size=13,face="bold"),
      plot.tag.position=c(0.01,0.98),
      legend.position="none",
      plot.margin=margin(0,0,bottom_space,0)
    )
  
  if(show_x_tick){
    p<-p+
      theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_line(color="black",linewidth=0.5)
      )
  }else{
    p<-p+
      theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
      )
  }
  
  p
}

P_north_backcolor<-plot_slope_ts(
  slope_data,
  "Northern",
  hemi_cols["Northern"],
  "(a)",
  show_top_bar=TRUE,
  show_x_tick=TRUE,
  bottom_space=8
)

P_south_backcolor<-plot_slope_ts(
  slope_data,
  "Southern",
  hemi_cols["Southern"],
  "(b)",
  show_top_bar=FALSE,
  show_x_tick=FALSE,
  bottom_space=0
)

P_geo<-ggplot(data.frame(x=c(x_max_val,0),y=c(0,0)),aes(x=x,y=y))+
  geom_blank()+
  scale_x_reverse(
    name="Time (Ma)",
    limits=c(x_max_val,0),
    breaks=seq(500,0,-50),
    expand=c(0,0)
  )+
  scale_y_continuous(limits=c(0,0.01),expand=c(0,0))+
  # coord_geo(
  #   xlim=c(x_max_val,0),
  #   pos="bottom",
  #   dat=list("periods","epochs"),
  #   height=unit(1.45,"lines"),
  #   expand=FALSE
  # )+
  coord_geo(xlim =c(x_max_val,0),
            pos = as.list(rep("bottom", 2)),
            # dat = list("periods", "era"),
            dat = list("periods", era_uk),
            height = list(unit(1.35, "lines"), unit(1.35, "lines")),
            #fill="darkgrey",
            lab_color="black",
            rot = list(0, 0), 
            # size = list(4, 4), 
            abbrv =list(TRUE,FALSE)) +
  labs(x="Time (Ma)",y=NULL,tag=NULL)+
  theme_minimal()+
  theme(
    panel.grid=element_blank(),
    panel.border=element_blank(),
    panel.background=element_blank(),
    axis.title.x=element_text(size=14,margin=margin(t=2)),
    axis.text.x=element_text(size=14),
    axis.ticks.x=element_line(color="black",linewidth=0.5),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.tag=element_blank(),
    legend.position="none",
    plot.margin=margin(0,0,10,0)
  )
climate_bar<-ggplot(
  data.frame(climate_state=factor(climate_levels,levels=climate_levels)),
  aes(x=climate_state,y=1,fill=climate_state)
)+
  geom_tile(width=0.8,height=0.8)+
  scale_fill_manual(values=climate_colors,name="Climate state",drop=FALSE)+
  guides(fill=guide_legend(
    nrow=1,
    override.aes=list(alpha=1,colour=NA)
  ))+
  theme_void()+
  theme(
    legend.position="top",
    legend.title=element_text(size=12,face="bold"),
    legend.text=element_text(size=12),
    legend.margin=margin(4,0,4,0),
    legend.box.margin=margin(0,0,0,0),
    plot.margin=margin(0,0,0,0)
  )

climate_legend_grob<-cowplot::get_legend(climate_bar)

geo_climate_spacer<-cowplot::ggdraw()

aligned_plots<-cowplot::align_plots(
  P_north_backcolor,
  P_south_backcolor,
  P_geo,
  align="v",
  axis="lr"
)

slope_vTime_plot_backcolor<-cowplot::plot_grid(
  aligned_plots[[1]],
  aligned_plots[[2]],
  aligned_plots[[3]],
  geo_climate_spacer,
  climate_legend_grob,
  ncol=1,
  rel_heights=c(10,10,1.65,0.45,1.15)
) +
  theme(
  plot.margin=margin(10,10,10,10)
)

print(slope_vTime_plot_backcolor)

ts_jpg<-sprintf(
  "./figures/jpg/background color %skm %squota %s equal-area latitude bins per-cell balanced OLS time series.jpg",
  params$spacing,params$level,rich_params$n_lat_bins
)
ts_pdf<-sprintf(
  "./figures/pdf/background color %skm %squota %s equal-area latitude bins per-cell balanced OLS time series.pdf",
  params$spacing,params$level,rich_params$n_lat_bins
)
ggsave(
  filename=ts_jpg, plot=slope_vTime_plot_backcolor,
  width=8, height=7.4,  dpi=300, device="jpeg"
)
ggsave(
  filename=ts_pdf, plot=slope_vTime_plot_backcolor,
  width=8, height=7.4, dpi=300, device='pdf'
)
# -----------------------------------------------------------------------
# 6.Boxplot by climate state
# -----------------------------------------------------------------------
slope_data_filtered <- slope_data %>%
  filter(!is.na(slope_value),
         climate_state %in% climate_levels,
         label != "bad") %>%
  mutate(
    climate_state = factor(climate_state, levels = climate_levels),
    slope_type = factor(slope_type, levels = c("Northern", "Southern"))
  )
sample_counts <- slope_data_filtered %>%
  count(climate_state, name = "count") %>%
  mutate(
    climate_state = factor(climate_state, levels = climate_levels),
    state_with_n = paste0(climate_state, "\n(n=", count, ")")
  ) %>%
  arrange(climate_state)

slope_data_filtered <- slope_data_filtered %>%
  left_join(
    select(sample_counts, climate_state, state_with_n), by = "climate_state"
  ) %>%
  mutate(
    state_with_n = factor(state_with_n, levels = sample_counts$state_with_n))

y_min_box <- min(slope_data_filtered$slope_value, na.rm = TRUE)
y_max_box <- max(slope_data_filtered$slope_value, na.rm = TRUE)

boxplot <- ggplot(slope_data_filtered,
  aes(state_with_n, slope_value, fill = slope_type)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.1, ymax = 0.1,
    fill = "lightblue", alpha = 0.3) +
  geom_hline(
    yintercept = 0, color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_boxplot(
    outlier.shape = NA, position = position_dodge(0.75)) +
  geom_point(
    shape = 21,
    size = 1,
    alpha = 0.6,
    color = "black",
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.75
    ),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = hemi_cols, name = "Hemisphere") +
  labs(x = "Climate state", y = "Slope value") +
  coord_cartesian(clip = "off", xlim = c(1, 5), ylim = c(y_min_box, y_max_box)) +
  annotate(
    "text", x = 5.81, y = (y_max_box + 0.1) / 2,
    label = "Non-modern-type", size = 4.5, angle = 270) +
  annotate(
    "text", x = 5.81, y = (y_min_box - 0.1) / 2,
    label = "Modern-type", size = 4.5, angle = 270) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(10, 30, 10, 10)
  )

print(boxplot)

box_jpg<-sprintf(
  "./figures/jpg/%skm %squota %s equal-area latitude bins per-cell balanced OLS boxplot.jpg",
  params$spacing,params$level,rich_params$n_lat_bins
)
box_pdf<-sprintf(
  "./figures/pdf/%skm %squota %s equal-area latitude bins per-cell balanced OLS boxplot.pdf",
  params$spacing,params$level,rich_params$n_lat_bins
)

ggsave(box_jpg,boxplot,width=7,height=5,dpi=300)
ggsave(box_pdf,boxplot,width=7,height=5,dpi=300)

# -----------------------------------------------------------------------
# 7.Summary
# -----------------------------------------------------------------------

slope_summary<-slope_cli_df_filter%>%
  group_by(hemisphere)%>%
  summarise(
    good_count=n(),
    greater_than_0=sum(slope>0,na.rm=TRUE),
    less_than_0=sum(slope<0,na.rm=TRUE),
    equal_to_0=sum(slope==0,na.rm=TRUE),
    reverse_LDG=greater_than_0/good_count*100,
    normal_LDG=less_than_0/good_count*100,
    flat_LDG=equal_to_0/good_count*100,
    median_slope=median(slope,na.rm=TRUE),
    mean_slope=mean(slope,na.rm=TRUE),
    sd_slope=sd(slope,na.rm=TRUE),
    .groups="drop"
  )

summary_path<-sprintf(
  "./results/%skm %squota %s equal-area latitude bins per-cell balanced OLS slope summary.csv",
  params$spacing,params$level,rich_params$n_lat_bins
)
write.csv(slope_summary,summary_path,row.names=FALSE)

print(slope_summary)
print(slope_cli_path)
print(wil_path)
print(summary_path)