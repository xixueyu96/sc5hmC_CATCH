library(magrittr)
library(dplyr)
library(ggplot2)

# plot_dt <- readRDS("data/processed/DhMR.sc_5k.rds")
#
# mc_df <- read.table(
#   "data/source/WGBS_GSR.5000_1X.merged.bed",
#   header = T, stringsAsFactors = F
# )

plot_5mc_5hmC_ratio <- function(stage_pair){

  cmp_list <- colnames(mc_df)[-c(1:3)]
  names(cmp_list) <- c("ZY", "2C", "4C", "8C", "Morula", "Blast")

  pre_stage <- cmp_list[strsplit(stage_pair, "->")[[1]][1]]
  pro_stage <- cmp_list[strsplit(stage_pair, "->")[[1]][2]]

  in_fn <- file.path(
    "data/processed", paste0(
      gsub("->", "to", stage_pair),
      ".DhMR_5mC_ratio.WGBS_GSR.rds"
    )
  )

  tmp_df <- readRDS(in_fn)

  tmp_df %<>%
    filter(get(pre_stage)>0 & meth_g1>0 & meth_g2 > 0 & get(pro_stage) >0) %>%
    mutate(
      ratio_5hmC = log2((meth_g2) / (meth_g1) ), ## 5hmC
      ratio_5mC = log2((get(pro_stage)) / (get(pre_stage)) )) %>% ## 5mC
    mutate(
      demeth = case_when(
      ratio_5hmC > -1 & ratio_5mC < -1 ~ "active_de",
      TRUE ~ "other")
    )

  gc()

  # library(ggpubr)

  fill_pal <- c("stable"="grey", "decreasing"="#74b1d8", "increasing"="#df6c77")
  fill_pal <- c("other"="grey", "active_de"="#df6c77")

  plot1 <- ggplot(tmp_df, aes(x = ratio_5hmC, y = ratio_5mC, color = category)) +
    # geom_point(aes(color = demeth), size = 0.1, shape=21) +
    geom_point(aes(color = category), size = 0.1, shape=21) +
    # geom_density2d(aes(color=demeth), size=1) +
    geom_density2d(aes(color=category), size=1) +
    # geom_rug() +
    # geom_hline(yintercept = log2(1+1), linetype="longdash") +
    geom_hline(yintercept = log2(0.5), linetype="longdash") +
    # geom_vline(xintercept = c(log2(1.5), log2(3)), linetype="longdash") +
    geom_vline(xintercept = log2(0.5), linetype="longdash") +
    labs(y="logFC of 5mC", x="logFC of 5hmC") +
    scale_color_manual(values = fill_pal) +
    ylim(c(-2.5,2.5)) +
    labs(title = stage_pair) +
    theme_bw(base_size = 15) + theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = .5))

  p2 <- ggExtra::ggMarginal(plot1, groupColour = TRUE, groupFill = TRUE)
  print(p2)

  plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]

  ggsave(
    plot = p2, width = 5, height = 5, units = "in", dpi=320,
    filename = paste(
      "viz/DhMR_5mC_ratio.WGBS_GSR",
      plot_time,
      gsub("->", "to", stage_pair),
      "pdf", sep = ".")
  )

  cor_res <- cor.test(tmp_df$ratio_5hmC, tmp_df$ratio_5mC)
  print(cor_res)

}

lapply(
  c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),
  plot_5mc_5hmC_ratio
)
