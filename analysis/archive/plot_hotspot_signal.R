
plot_hotspot_hmC <- function(hotspot_rds, y_axis_title, plot_name){

  require(dplyr)
  require(ggplot2)
  require(ggridges)

  merged_gr <- readRDS(hotspot_rds)

  region_pal <- c("hotspot"="#7ca9d5", "matched"="#e7df68", "other"="#e5545b")

  merged_long <- merged_gr %>%
    as.data.table() %>%
    dplyr::select(`X2C`, `X4C`, `X8C`, `Morula`, `Blast`, `region`) %>%
    reshape2::melt(id.var="region") %>%
    filter(!is.na(value))

  merged_long$region <- factor(
    merged_long$region,
    levels = c("other", "matched", "hotspot")
  )

  p <- merged_long %>%
    ggplot(aes(x=region, y=value, color=region))+
    geom_jitter(aes(group=region),shape=21, width = .3) +
    geom_violin(scale = "width", fill=NA, color="black") +
    facet_wrap(.~variable, nrow = 2, scales = "free") +
    scale_x_discrete(guide = "axis_nested") +
    scale_color_manual(values = region_pal) +
    labs(y=y_axis_title) +
    egg::theme_presentation() +
    theme(aspect.ratio = 1,
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

  ggsave(
    plot = p,
    filename = plot_name,
    width = 12,
    height = 6
  )

  # return(p)

}

plot_hotspot_hmC(
  "data/processed/hotspot.5hmC.1X50.240218.rds",
  y_axis_title = "5hmCpG level",
  plot_name = "viz/hotspot.5hmC.1X50.240218.pdf"
)
plot_hotspot_hmC(
  "data/processed/hotspot.5hmC_5mC_nor.1X50.240218.rds",
  y_axis_title = "5hmCpG/5mCpG level",
  plot_name = "viz/hotspot.5hmC_5mC_nor.1X50.240218.pdf"
)
