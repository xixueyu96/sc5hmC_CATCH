library(dplyr)

plot_dt <- readRDS("data/processed/DhMR.sc_5k.rds")

mc_df <- read.table(
  "data/source/WGBS_GSR.5000_1X.merged.bed",
  header = T, stringsAsFactors = F
)

calc_5mc_DhMR <- function(stage_pair){

  cmp_list <- colnames(mc_df)[-c(1:3)]
  names(cmp_list) <- c("ZY", "2C", "4C", "8C", "Morula", "Blast")

  pre_stage <- cmp_list[strsplit(stage_pair, "->")[[1]][1]]
  pro_stage <- cmp_list[strsplit(stage_pair, "->")[[1]][2]]

  # stage_pair <- "Morula->Blast"

  tmp_df <- plot_dt %>%
    filter(cmp==stage_pair) %>%
    inner_join(
      # mc_df %>% mutate(End=Pos+5e3-1),
      mc_df %>%
        mutate(chr=paste0("chr", chrom)) %>%
        dplyr::select(
          "chr", "end",
          all_of(as.character(pre_stage)),
          all_of(as.character(pro_stage))),
      by = c("chr","end")) %>%
    rename(
      !!"meth_mc1" := !!sym(pre_stage),
      !!"meth_mc2" := !!sym(pro_stage)
    )

  out_fn <- file.path(
    "data/processed", paste0(
      gsub("->", "to", stage_pair),
      ".DhMR_5mC_ratio.WGBS_GSR.rds"
    )
  )

  saveRDS(tmp_df, out_fn)

}

lapply(
  c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),
  calc_5mc_DhMR
)
