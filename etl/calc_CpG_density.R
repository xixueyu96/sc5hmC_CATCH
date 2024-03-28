##--------------coordinate of CpG sites----------------------

library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)

chrs <- names(Mmusculus)[1:22]
cgs <-
  lapply(chrs, function(x)
    start(matchPattern("CG", Mmusculus[[x]])))
gcs <-
  lapply(chrs, function(x)
    start(matchPattern("GC", Mmusculus[[x]])))
cpgr <-
  do.call(c, lapply(1:22, function(x)
    GRanges(
      names(Mmusculus)[x], IRanges(cgs[[x]], width = 2)
    )))
gpcr <-
  do.call(c, lapply(1:22, function(x)
    GRanges(
      names(Mmusculus)[x], IRanges(gcs[[x]], width = 2)
    )))
end(cpgr) <- start(cpgr)
start(gpcr) <- end(gpcr)
all_cpg <- union(cpgr, gpcr)
all_cpg$CpG <- 1

all_cpg

saveRDS(all_cpg, file = "data/handmade/mm10.CpG.rds")

gc()

##-------------- coordinate of Ns ----------------------
chrs <- names(Mmusculus)[1:22]
ns <-
  lapply(chrs, function(x)
    start(matchPattern("N", Mmusculus[[x]])))
n_gr <-
  do.call(c, lapply(1:21, function(x)
  {
    GRanges(names(Mmusculus)[x], IRanges(ns[[x]], width = 2))
  }))
n_gr <- reduce(n_gr)

n_gr

saveRDS(n_gr, "data/handmade/mm10.Ns.rds")

gc()

##--------------CpG density in tiles----------------------
calc_cpg_density_tile <- function(window_size){
  require(GenomicRanges)
  require(circlize)
  require(EnrichedHeatmap)

  # setwd("/mnt/f/Project/sc5hmC/data")
  all_cpg <- readRDS("data/handmade/mm10.CpG.rds")
  n_gr <- readRDS("data/handmade/mm10.Ns.rds")

  # window_size <- 1e4

  chr_df <-  read.chromInfo(species="mm10")$df
  chr_gr <- GRanges(
    seqnames = chr_df[, 1],
    ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3])
  )
  chr_window <- makeWindows(chr_gr, w = window_size)

  ## filter Ns
  hits <- overlapsAny(chr_window, n_gr)
  chr_noN_windows <- chr_window[!hits]
  # chr_noN_windows

  ## from jzh @ biubiu-1s
  seqlevels(all_cpg)
  # seqlevels(chr_window)
  seqlevels(chr_noN_windows)

  all_cpg <- sort(all_cpg)
  # chr_window <- sort(chr_window)
  chr_noN_windows <- sort(chr_noN_windows)

  hits <- findOverlaps(all_cpg, chr_noN_windows)
  hits_rle <- rle(sort(hits@to))
  chr_noN_windows$cpg_density <- 0
  chr_noN_windows$cpg_density[hits_rle$values] <- hits_rle$lengths/window_size

  gc()
  # chr_noN_windows
  return(chr_noN_windows)
}

chr_noN_windows <- calc_cpg_density_tile(1e3)

hist_df <- hist(chr_noN_windows$cpg_density, breaks = 100)

## save binned cpg_density
saveRDS(chr_noN_windows, file = "data/processed/mm10.CpG_density.win1k.noN.rds")

## --------------CpG density around detected site----------------------

calc_cpg_density_center <- function(hmC_bed, out_fn){
  require(GenomicRanges)
  require(data.table)
  require(dplyr)

  window_size <- 1000

  # tet_cpg_dt <- fread("cut -f 1,3-5 data/source/TET3_cKO_site/TET3_KO_2C.1X50.aim.bed") %>%
  #   mutate(V5=V2) %>% setNames(c("seqnames", "start", "met", "total", "end"))
  # cmd <- paste0("cut -f 1,3-5 ", hmC_bed)
  all_cpg <- readRDS("data/handmade/mm10.CpG.rds")

  tet_cpg_dt <- fread(hmC_bed) %>% mutate(V2=V3) %>%
    setNames(c("seqnames", "start", "end", "met", "unmet")) %>% filter(met > 0)
  tet_cpg_gr <- makeGRangesFromDataFrame(tet_cpg_dt, keep.extra.columns = T)
  tet_cpg_gr$met_ratio <- tet_cpg_gr$met/(tet_cpg_gr$unmet + tet_cpg_gr$met)
  # n_gr <- readRDS("data/handmade/mm10.Ns.rds")

  tet_window <- resize(tet_cpg_gr, width = window_size, fix = "center")

  ## from jzh @ biubiu-1s
  seqlevels(tet_cpg_gr)
  # seqlevels(chr_window)
  seqlevels(tet_window)

  all_cpg <- sort(all_cpg)
  # chr_window <- sort(chr_window)
  tet_window <- sort(tet_window)

  hits <- findOverlaps(all_cpg, tet_window)

  gc()

  hits_rle <- rle(hits@to)
  tet_window$cpg_density <- 0
  tet_window$cpg_density[hits_rle$values] <- hits_rle$lengths/window_size

  gc()

  saveRDS(tet_window, out_fn)
}

calc_cpg_density_center(
  "data/source/TET3_cKO_site/TET3_KO_2C.1X50.aim.bed",
  "data/processed/TET3_KO_2C.1kb_cpg_density.rds"
)
calc_cpg_density_center(
  "data/source/TET3_cKO_site/TET3_KO_EZ.1X50.aim.bed",
  "data/processed/TET3_KO_EZ.1kb_cpg_density.rds"
)
calc_cpg_density_center(
  "data/source/TET3_cKO_site/TET3_KO_LZ.1X50.aim.bed",
  "data/processed/TET3_KO_LZ.1kb_cpg_density.rds"
)

