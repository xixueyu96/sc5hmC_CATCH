library(data.table)
library(dplyr)

##-----------------------5hmC-----------------------

fn <- list.files(
  "data/source/hmC_tile/",
  pattern = ".1X50.aim.1000_1X.bed",
  full.names = T
)

df_list <- lapply(
  fn, function(x){

    df <- fread(x, stringsAsFactors = F)

    sp <- strsplit(basename(x), "[.]")[[1]][1]

    colnames(df) <- c("seqnames", "start", "end", "cov", "meth")

    df <- df %>%
      mutate(ratio=meth/cov) %>%
      dplyr::select(
        seqnames, start, end, ratio)

    colnames(df) <- c("seqnames", "start", "end", sp)

    return(df)
  }
)

mtx <- df_list %>%
  purrr::reduce(
    function(x, y){
      base::merge(x, y, by = c("seqnames", "start", "end"), all = T)
    }) %>%
  dplyr::select(
    `seqnames`, `start`, `end`,
    `Oocyte`, `Sperm`, early_zygote, late_zygote,
    `2C`, `4C`, `8C`, `Morula`, `Blast`
  )

rm(df_list)

saveRDS(mtx, file = "data/processed/mouse_embryo.5hmC.1X50.aim.1kb.rds")

##----------------------- 5mC -----------------------

fn <- list.files(
  "data/source/hmC_mC_nor_tile/",
  pattern = ".hmc_mc.1k_tile.bedGraph",
  full.names = T
)

df_list <- lapply(
  fn, function(x){

    df <- fread(x, stringsAsFactors = F)

    sp <- strsplit(basename(x), "[.]")[[1]][1]
    sp <- strsplit(sp, "_")[[1]][2]

    colnames(df) <- c("seqnames", "start", "end", "ratio")

    df <- df %>%
      mutate(start = start + 1) %>%
      filter(!is.na(ratio))

    colnames(df) <- c("seqnames", "start", "end", sp)

    return(df)
  }
)

mtx <- df_list %>%
  purrr::reduce(
    function(x, y){
      base::merge(x, y, by = c("seqnames", "start", "end"), all = T)
    }) %>%
  dplyr::select(
    `seqnames`, `start`, `end`,
    `Oocyte`, `Sperm`, ZY,
    `2C`, `4C`, `8C`, `Morula`, `Blast`
  )

rm(df_list)

saveRDS(mtx, file = "data/processed/mouse_embryo.5hmC_5mC_nor.1X50.aim.1kb.rds")
