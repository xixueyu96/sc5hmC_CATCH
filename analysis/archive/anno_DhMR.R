##---------------------active DhMR (rGREAT)---------------------

library(rGREAT)
library(GenomicRanges)
library(dplyr)

lapply(

  c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),

  function(x){
    # stage_pair <- "2C->4C"

    in_fn <- file.path(
      "data/processed", paste0(
        gsub("->", "to", x),
        ".DhMR_5mC_ratio.WGBS_GSR.rds"
      )
    )

    tmp_df <- readRDS(in_fn)

    sub_gr <- tmp_df %>%
      filter(meth_g2 > 0.5 * meth_g1) %>%
      dplyr::select(chr, start, end) %>%
      makeGRangesFromDataFrame()

    job <- submitGreatJob(
      sub_gr,species="mm10",
      rule = "basalPlusExt",
      request_interval = 50
    )

    # res_great <- plotRegionGeneAssociationGraphs(job)
    tb_great <- getEnrichmentTables(job)
    # View(tb_great$`GO Biological Process`)
    ## save result
    out_fn <- file.path(
      "data/processed",
      paste0(gsub("->", "to", x),".DhMR_GREAT.tsv"
      ))

    write.table(
      tb_great$`GO Biological Process`, out_fn,
      col.names = T, row.names = F, sep = "\t", quote = F
    )
  }
)


##---------------------active DhMR (LOLA)---------------------
library(circlize)
library(EnrichedHeatmap)
library(LOLA)

## generate control dataset
chr_df <-  read.chromInfo(species="mm10")$df
chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_window <- makeWindows(chr_gr, w = 5e3)

## reference dataset
regionDB <- loadRegionDB(
  "/mnt/f/packages/LOLACore/mm10",
  collections = "codex"
)

lapply(

  c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),

  function(x){

    # stage_pair <- "2C->4C"
    in_fn <- file.path(
      "data/processed", paste0(
        gsub("->", "to", x),
        ".DhMR_5mC_ratio.WGBS_GSR.rds"
      )
    )

    tmp_df <- readRDS(in_fn)

    sub_gr <- tmp_df %>%
      filter(meth_g2 > 0.5 * meth_g1) %>%
      dplyr::select(chr, start, end) %>%
      makeGRangesFromDataFrame()

    userSets <- GRangesList(sub_gr)
    names(userSets) <- x

    ## calc enrichment
    locResults <- runLOLA(
      userSets, chr_window, regionDB,
      minOverlap = 1, direction = "enrichment"
    )

    ## save result
    out_fn <- file.path(
      "data/processed",
      paste0(gsub("->", "to", x),".DhMR_LOLA.tsv"
      ))

    write.table(
      locResults, out_fn,
      col.names = T, row.names = F, sep = "\t", quote = F
    )
  }
)
