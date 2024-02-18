library(LOLA)
library(dplyr)
library(ggplot2)

## load input
control_gr <- readRDS("data/processed/hotspot.cpg_density_matched.rds")
hotspot_gr <- readBed("data/source/1111_hmC_ratio_1X50.merge.1k.v5s4.hotspot.bed")

run_lola <- function(x){

  regionDB <- loadRegionDB(
    "/mnt/f/packages/LOLACore/mm10/",
    collections = x
  )

  userSets <- GRangesList(hotspot_gr)
  names(userSets) <- "hotspot"

  locResults <- runLOLA(
    userSets,
    c(hotspot_gr, control_gr),
    regionDB,
    minOverlap = 1,
    direction = "enrichment"
  )

  plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
  fn_1 <- paste0("data/processed/hotspot.lola_", x, ".", plot_time, ".rds")
  saveRDS(locResults, fn_1)

  p <- locResults %>%
    mutate(label=ifelse(qValue < 0.05, description, "")) %>%
    ggplot(aes(x=oddsRatio, y=-log10(qValue)))+
    geom_point(size=2, shape=21, fill="#f79c56")+
    ggrepel::geom_text_repel(
      aes(label=label),box.padding = 0.5, max.overlaps =5)+
    theme_classic(base_size = 20) +
    theme(aspect.ratio = 1)

  fn_2 <- paste0("viz/hotspot.lola_", x, ".", plot_time, ".pdf")

  ggsave(
    plot = p,
    filename = fn_2,
    width = 10,
    height = 10,
    units = "cm"
  )
}

run_lola("rmsk_gene")
run_lola("rmsk_LINE")
run_lola("rmsk_LTR")
