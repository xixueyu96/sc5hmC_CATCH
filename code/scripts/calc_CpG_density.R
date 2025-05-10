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

saveRDS(all_cpg, file = "data/metadata/mm10.CpG.rds")

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

saveRDS(n_gr, "data/metadata/mm10.Ns.rds")

gc()

##--------------CpG density in tiles----------------------
calc_cpg_density_tile <- function(window_size){
  require(GenomicRanges)
  require(circlize)
  require(EnrichedHeatmap)

  # setwd("/mnt/f/Project/sc5hmC/data")
  all_cpg <- readRDS("data/metadata/mm10.CpG.rds")
  n_gr <- readRDS("data/metadata/mm10.Ns.rds")

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
saveRDS(chr_noN_windows, file = "data/metadata/mm10.CpG_density.win1k.noN.rds")
