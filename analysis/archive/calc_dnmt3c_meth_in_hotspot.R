library(GenomicRanges)
library(EnrichedHeatmap)
library(rtracklayer)
library(data.table)

region_pal <- c("hotspot"="#7ca9d5", "matched"="#e7df68", "other"="#e5545b")

merged_gr <- readRDS("data/processed/hotspot.5hmC.1X50.240218.rds")

germ_het <- import("data/public/GSE94688_RPM_log2_500bp_tiles_ChIP.txt.gz")
germ_het

mat1 <-  normalizeToMatrix(
  germ_het, merged_gr, value_column = "score",
  # background = NA, smooth = FALSE,
  extend = 3000, mean_mode = "absolute", w = 50
)
mat1
# rm(pol2_E2C)

top_anno <- HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = region_pal)))
meth_col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

ht1 <- EnrichedHeatmap(
  mat1, col = meth_col_fun,
  name = "gem_meth",
  top_annotation = top_anno,
  row_split = merged_gr$region,
  column_title = "meth in hotspot",
  # width=unit(4,"cm"),
  use_raster=T,
  raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)
ht1
