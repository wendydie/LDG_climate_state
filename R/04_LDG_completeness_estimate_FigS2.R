# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 04_LDG_completeness_estimate_FigS2.R
# Last updated: 2026-08-18
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(palaeoverse); library(grid); library(cowplot)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

# -----------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------
rich_df <- read.csv(sprintf(
  "./results/LDG/%s_cell_%s_richness.csv",
  params$spacing, params$level
))
time_bins <- readRDS("./data/time_bins.RDS")
lat_bins <- palaeoverse::lat_bins_area(n=rich_params$n_lat_bins) %>% arrange(min)

rich_df <- rich_df %>%
  mutate(stage=time_bins$interval_name[match(bin_midpoint,time_bins$mid_ma)]) %>%
  filter(bin_midpoint<=486.85,!is.na(stage)) %>%
  mutate(
    completeness=ifelse(nT>=5 & t<=2*nT,"Complete","Incomplete"),
    hemisphere=ifelse(cell_lat>=0,"Northern","Southern"),
    bin_index=findInterval(cell_lat,c(lat_bins$min,Inf)),
    bin=lat_bins$bin[bin_index]
  ) %>%
  left_join(lat_bins %>% select(bin,lat_bin_mid=mid),by="bin") %>%
  mutate(abs_lat_bin_mid=abs(lat_bin_mid))

# -----------------------------------------------------------------------
# 2. QC
# -----------------------------------------------------------------------
rich_complete <- rich_df %>%
  filter(completeness=="Complete") %>%
  mutate(
    lat_zone=case_when(
      abs_lat_bin_mid<30~"tropical",
      abs_lat_bin_mid<60~"temperate",
      abs_lat_bin_mid<=90~"polar",
      TRUE~NA_character_
    )
  )

adjacent_df <- has_adjacent_bins(rich_complete,lat_bins) %>%
  distinct(bin_midpoint,hemisphere,label) %>%
  transmute(bin_midpoint,hemisphere,has_adjacent_tt=label=="good")

qc <- rich_complete %>%
  group_by(bin_midpoint,hemisphere) %>%
  summarise(
    has_tropical=any(lat_zone=="tropical",na.rm=TRUE),
    has_temperate=any(lat_zone=="temperate",na.rm=TRUE),
    .groups="drop"
  ) %>%
  left_join(adjacent_df,by=c("bin_midpoint","hemisphere")) %>%
  mutate(
    has_adjacent_tt=coalesce(has_adjacent_tt,FALSE),
    label=ifelse(has_tropical & has_temperate & has_adjacent_tt,"good","bad")
  )

rich_df2 <- rich_df %>%
  left_join(qc %>% select(bin_midpoint,hemisphere,label),
            by=c("bin_midpoint","hemisphere")) %>%
  mutate(
    label=coalesce(label,"bad"),
    completeness=factor(completeness,levels=c("Incomplete","Complete"))
  )

# -----------------------------------------------------------------------
# 3. Era assignment
# -----------------------------------------------------------------------
stage_order <- rich_df2 %>%
  distinct(bin_midpoint,stage) %>%
  arrange(desc(bin_midpoint)) %>%
  mutate(id=row_number())

req <- c("Tremadocian","Changhsingian","Induan","Maastrichtian","Danian")
miss <- setdiff(req,stage_order$stage)
if(length(miss)>0) stop("Missing stage(s): ",paste(miss,collapse=", "))

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

rich_df2 <- rich_df2 %>%
  left_join(stage_order %>% select(bin_midpoint,era),by="bin_midpoint")

# -----------------------------------------------------------------------
# 4. Plot function
# -----------------------------------------------------------------------
make_era_plot <- function(era_name,show_x=FALSE){
  
  lev <- stage_order %>%
    filter(era==era_name) %>%
    arrange(desc(bin_midpoint)) %>%
    pull(stage)
  
  df <- rich_df2 %>%
    filter(era==era_name) %>%
    mutate(stage=factor(stage,levels=lev))
  
  bad_bg <- df %>%
    filter(label=="bad") %>%
    distinct(stage,hemisphere) %>%
    mutate(
      xmin=ifelse(hemisphere=="Northern",0,-90),
      xmax=ifelse(hemisphere=="Northern",90,0),
      text_x=ifelse(hemisphere=="Northern",45,-45)
    )
  
  ggplot(df,aes(cell_lat,fill=completeness))+
    geom_rect(
      data=bad_bg,
      aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
      fill="gray90",alpha=.3,inherit.aes=FALSE
    )+
    geom_histogram(
      binwidth=5,position="stack",
      color="black",linewidth=.25
    )+
    geom_vline(xintercept=0,color="red",linewidth=.7)+
    geom_text(
      data=bad_bg,
      aes(x=text_x,y=Inf,label="Bad"),
      color="red",size=2.5,fontface="bold",vjust=1.25,
      inherit.aes=FALSE
    )+
    scale_fill_manual(
      name="Completeness",
      values=c(Incomplete="#F0E442",Complete="#009E73")
    )+
    scale_x_continuous(
      breaks=seq(-90,90,30),
      expand=c(0,0)
    )+
    coord_cartesian(xlim=c(-90,90))+
    scale_y_continuous(
      breaks=function(y){
        m <- ceiling(max(y,na.rm=TRUE)/10)*10
        unique(round(seq(0,m,length.out=5)))
      }
    )+
    facet_wrap(~stage,ncol=6,scales="free_y",drop=FALSE)+
    labs(
      x=if(show_x) "Palaeolatitude" else NULL,
      y=NULL
    )+
    theme_minimal()+
    theme(
      strip.text=element_text(
        size=6.8,face="bold",
        margin=margin(0,0,0,0)
      ),
      strip.placement="inside",
      
      panel.spacing.x=unit(.06,"lines"),
      panel.spacing.y=unit(0,"lines"),
      
      panel.border=element_rect(
        color="black",fill=NA,linewidth=.65
      ),
      panel.grid=element_blank(),
      panel.background=element_rect(fill="white"),
      
      axis.ticks=element_line(color="black",linewidth=.4),
      axis.ticks.length=unit(.035,"cm"),
      axis.text=element_text(size=6.8,color="black"),
      axis.title.x=element_text(size=9.5,color="black",margin=margin(t=2)),
      
      legend.position="none",
      
      plot.background=element_rect(
        color="grey40",
        fill=NA,
        linewidth=.7,
        linetype="dashed"
      ),
      
      # facet 与虚线框距离
      plot.margin=margin(1.5,1.5,1.5,1.5)
    )
}

# -----------------------------------------------------------------------
# 5. Era labels
# -----------------------------------------------------------------------
era_label <- function(x){
  ggplot()+
    annotate(
      "text",x=.5,y=.5,label=x,
      angle=270,fontface="bold",size=3.2
    )+
    xlim(0,1)+ylim(0,1)+
    theme_void()+
    theme(plot.margin=margin(0,0,0,0))
}

# -----------------------------------------------------------------------
# 6. Three Era blocks
# -----------------------------------------------------------------------
p_pal <- make_era_plot("Palaeozoic",FALSE)
p_mes <- make_era_plot("Mesozoic",FALSE)
p_cen <- make_era_plot("Cenozoic",TRUE)

row_pal <- p_pal + era_label("Palaeozoic") +
  plot_layout(widths=c(1,.018))

row_mes <- p_mes + era_label("Mesozoic") +
  plot_layout(widths=c(1,.018))

row_cen <- p_cen + era_label("Cenozoic") +
  plot_layout(widths=c(1,.018))

# -----------------------------------------------------------------------
# 7. Combine Era blocks
# -----------------------------------------------------------------------
n_pal <- ceiling(sum(stage_order$era=="Palaeozoic")/6)
n_mes <- ceiling(sum(stage_order$era=="Mesozoic")/6)
n_cen <- ceiling(sum(stage_order$era=="Cenozoic")/6)

p_body <-
  row_pal /
  plot_spacer() /
  row_mes /
  plot_spacer() /
  row_cen +
  plot_layout(
    heights=c(
      n_pal,.008,
      n_mes,.008,
      n_cen
    )
  )

# -----------------------------------------------------------------------
# 8. External legend
# -----------------------------------------------------------------------
legend_plot <- ggplot(
  rich_df2,
  aes(cell_lat,fill=completeness)
)+
  geom_histogram(binwidth=5)+
  scale_fill_manual(
    name="Completeness",
    values=c(
      Incomplete="#F0E442",
      Complete="#009E73"
    )
  )+
  guides(
    fill=guide_legend(
      nrow=1,
      title.position="left"
    )
  )+
  theme_void()+
  theme(
    legend.position="bottom",
    legend.direction="horizontal",
    legend.title=element_text(size=8.5),
    legend.text=element_text(size=7.5),
    legend.margin=margin(-2,0,-2,0),
    legend.box.margin=margin(0,0,0,0),
    legend.spacing.x=unit(2,"pt"),
    legend.key.height=unit(5,"pt"),
    legend.key.width=unit(12,"pt")
  )

legend_grob <- cowplot::get_legend(legend_plot)

# -----------------------------------------------------------------------
# 9. Final compact layout
# -----------------------------------------------------------------------

body_grob <- patchwork::patchworkGrob(p_body)

main_with_y <- cowplot::plot_grid(
  grid::textGrob(
    "Cell number",
    rot=90,
    gp=grid::gpar(fontsize=10.5)
  ),
  body_grob,
  nrow=1,
  rel_widths=c(.018,1)
)

p_final <- cowplot::plot_grid(
  main_with_y,
  legend_grob,
  ncol=1,
  rel_heights=c(1,.025)
)

# -----------------------------------------------------------------------
# 10. Save
# -----------------------------------------------------------------------
print(p_final)

fig_path <- sprintf(
  "./figures/FigS2_%s_km_%s_quota_%s.jpg",
  params$spacing,
  params$level,
  rich_params$n_lat_bins
)

ggsave(
  filename=fig_path,
  plot=p_final,
  width=8,
  height=8.5,
  dpi=300,
  bg="white"
)

cat("\nSaved: ",fig_path,"\n",sep="")