library(dplyr)
library(data.table)
library(GenomicRanges)

control_gr <- readRDS("data/processed/hotspot.cpg_density_matched.rds")
hotspot_gr <- LOLA::readBed("data/source/1111_hmC_ratio_1X50.merge.1k.v5s4.hotspot.bed")

calc_hotspot_signal <- function(fn_in, fn_out){

  fn_in <- paste0("data/processed/", fn_in)
  signal_tile <- fread(fn_in)

  merged_gr <- makeGRangesFromDataFrame(signal_tile, keep.extra.columns = T)

  seqlevels(merged_gr)
  merged_gr <- sort(merged_gr)

  hotspot_hit <- findOverlaps(hotspot_gr, merged_gr)
  control_hit <- findOverlaps(control_gr, merged_gr)

  merged_gr$region <- "other"
  merged_gr$region[hotspot_hit@to] <- "hotspot"
  merged_gr$region[control_hit@to] <- "matched"

  merged_other_gr <- merged_gr %>% .[.$region=="other"]
  merged_other_gr <- merged_other_gr[sample(length(merged_other_gr), size = length(hotspot_gr))]
  merged_gr <- c(merged_other_gr, merged_gr %>% .[.$region!="other"])
  merged_gr$region <- factor(merged_gr$region, levels = c("other", "matched", "hotspot"))
  merged_gr <- sort(merged_gr)

  saveRDS(merged_gr, file = paste0("data/processed/", fn_out))

  write.table(
    merged_gr[,"region"] %>% as.data.frame,
    file = paste0("data/processed/",gsub("rds", "bed", fn_out)),
    col.names = T, row.names = F, sep = "\t", quote = F
  )
}

calc_hotspot_signal(
  "mouse_embryo.5hmC.1X50.aim.1kb.230104.bed",
  "hotspot.5hmC.1X50.240218.rds"
)

calc_hotspot_signal(
  "mouse_embryo.5hmC_5mC_nor.1X50.aim.1kb.230104.bed",
  "hotspot.5hmC_5mC_nor.1X50.240218.rds"
)

