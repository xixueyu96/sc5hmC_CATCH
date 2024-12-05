
source("etl/merge_DHMR_lfc.R")
source("etl/plot_DHMR.frac.R")
source("etl/plot_RDL_RhML.R")

type_col_list <- c("stable"="#329cc8","changing"="#f6bb46")
category_col_list <- c("increasing"="#ee6300","decreasing"="#97ca5b")

cmp_list <- c("Sperm","LZY","2C","4C","8C","Morula")
names(cmp_list) <- c("EZY","2C","4C","8C","Morula", "Blast")
levels <- c(
  "Sperm->EZY","LZY->2C",
  "2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"
)

# plot_dt <- merge_DHMR_lfc("data/source/metilene_out/",cmp_list, "5k", min_lfc = 1)
# saveRDS(plot_dt, file = "data/processed/DhMR.sc_5k.rds")

plot_dt <- readRDS("data/processed/DhMR.sc_5k.rds")

p2 <- plot_DHMR.frac.stable(plot_dt, save_file=F, level=rev(levels))
p2 <- p2+theme(axis.text.y = element_blank())
p1 <- plot_DHMR.abs.stable(plot_dt, save_file=F, level=rev(levels))
p4 <- plot_DHMR.frac.changing(plot_dt, save_file=F, level=rev(levels))
p4 <- p4+theme(axis.text.y = element_blank())
p3 <- plot_DHMR.abs.changing(plot_dt, save_file=F, level=rev(levels))
p3 <- p3+theme(axis.text.y = element_blank())

pdf("viz/Fig4B.DhMR.5k.pdf", width = 10, height = 4)
egg::ggarrange(p1,p2,p3,p4,nrow = 1)
dev.off()

# ggsave("viz/Fig4B.DhMR.5k.pdf", width = 27.5, height = 24, units = "cm")


p5 <- plot_RDL.v2(plot_dt, save_file=F, level=rev(levels))
# p5 <- p5+theme(axis.text.y = element_blank())
p6 <- plot_RhML.v2(plot_dt, save_file=F, level=rev(levels))
p6 <- p6+theme(axis.text.y = element_blank())

pdf("viz/FigS4A.RhML_RDhML.5k.pdf", width = 5, height = 5)
egg::ggarrange(p5,p6,nrow = 1)
dev.off()


##---------------------plot 5hmC abundance in WT and TETcKO late zygote and 2-cell---------------------
LZto2C <- plot_dt %>% filter(cmp=="LZY->2C") %>%
  mutate(start = end - 5000) %>%
  select(chr, start, end, category)

options(scipen=999)
write.table(
  LZto2C, file = "data/processed/LZYto2C.DhmR.sc_5k.bed",
  quote = F, sep = "\t", col.names = F, row.names = F
)
options(scipen=2)
# bash calc_5hmC_in_DhmR.sh WT_2C.sample.xls

meta <- lapply(
  list.files("data/documentation", pattern = "sample.xls", full.names = T),
  function(x){
    df <- read.table(x) %>% setNames("Sample") %>%
      mutate(
        Group = gsub(".sample.xls", "", basename(x))
      )
  }
)

meta <- do.call("rbind", meta)

DhmR_5hmC_abundance <- read.table(
    "data/processed/WT_TETcKO.LZY_to_2C.DhmR.sc_5k.abundance.tsv",
    header = F, stringsAsFactors = F) %>%
  setNames(c("Sample", "category", "abundance")) %>%
  inner_join(meta, by = "Sample")

DhmR_5hmC_abundance$category <- factor(
  DhmR_5hmC_abundance$category,
  levels = c("stable", "decreasing", "increasing")
)

DhmR_5hmC_abundance$Group <- factor(
  DhmR_5hmC_abundance$Group,
  levels = c("WT_LZ", "TET3cKO_LZ", "WT_2C", "TET3cKO_2C")
)

df.summary <- DhmR_5hmC_abundance %>%
  group_by(category, Group) %>%
  summarise(
    # sd = sd(abundance, na.rm = TRUE),
    se = sd(abundance) / sqrt(length(abundance)),
    abundance = mean(abundance)
  )
df.summary

library(broom)
library(tidyr)
library(purrr)

df.test <- DhmR_5hmC_abundance %>%
  separate(
    col = "Group",
    sep = "_",
    into = c("Group", "Stage")) %>%
  nest(data = c(-category, -Stage)) %>%
  mutate(data = map(data, ~ wilcox.test(abundance ~ Group, data = .x)),
         data = map(data, tidy)) %>%
  unnest(data) %>%
  mutate(
    qvalue = p.adjust(p.value, "fdr")
  )
df.test

write.table(
  df.test, file = "data/processed/WT_TETcKO.LZY_to_2C.DhmR.sc_5k.abundance_test.tsv",
  col.names = T, row.names = F, sep = "\t", quote = F
)

df.summary %>%
  filter(category!="stable") %>%
  ggplot(aes(category, abundance, fill=Group)) +
  # geom_col(fill = "lightgray", color = "black") +
  geom_errorbar(
    aes(ymin = abundance-se, ymax = abundance+se),
    position = position_dodge(0.9), width=.5)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(y="5hmCpG levels", x="Dynamic DhMRs") +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 0.75)

ggsave(
  "viz/WT_TET3cKO.LZY_to_2C.DhmR.sc_5k.abundance.pdf",
  width = 8, height = 6
)

##-------------------- plot 5mC in DhMR------------------
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

  mc_DhMR_df <- readRDS(in_fn)

  mc_DhMR_df %<>%
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

  plot1 <- ggplot(mc_DhMR_df, aes(x = ratio_5hmC, y = ratio_5mC, color = category)) +
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

##---------------------flux of DhMR---------------------

# remotes::install_github("davidsjoberg/ggsankey")
# library(ggsankey)
library(ggplot2)
library(dplyr)
library(ggalluvial)

plot_dt <- readRDS("data/processed/DhMR.sc_5k.rds")

sankey_df <- plot_dt %>%
  # filter(category!="stable") %>%
  mutate(pos=paste0(chr, ":", start, "-", end)) %>%
  dplyr::select(pos, cmp, category) %>%
  tidyr::spread(key = pos, value = category) %>%
  tibble::column_to_rownames("cmp") %>% t

pal <- c("stable"="#329cc8","increasing"="#ee6300","decreasing"="#97ca5b", NA="grey")

sankey_long <- sankey_df %>%
  as.data.frame() %>%
  dplyr::select(`LZY->2C`,`2C->4C`, `4C->8C`, `8C->Morula`, `Morula->Blast`) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  group_by(`LZY->2C`,`2C->4C`, `4C->8C`, `8C->Morula`, `Morula->Blast`) %>%
  summarise(comb=n())

## convert NA to no_cover
sankey_long %>%
  filter(comb > 10) %>%
  ggplot(
    aes(
      # axis1 = `Sperm->EZY`,
      axis1 = `LZY->2C`,
      axis2 = `2C->4C`,
      axis3 = `4C->8C`,
      axis4 = `8C->Morula`,
      axis5 = `Morula->Blast`,
      y = comb)) +
  geom_flow(aes(fill=after_stat(stratum))) +
  geom_stratum(aes(fill=after_stat(stratum)), color="white",width = 1/4) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5),
    labels = c("LZY->2C", "2C->4C", "4C->8C", "8C->Morula", "Morula->Blast")) +
  # geom_text(stat = "stratum",
  #           aes(label = after_stat(stratum))) +
  theme_alluvial(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "viz/Fig4B.DhMR_flux.5k.pdf",
  width = 8, height = 4
)

plot_dt %>%
  filter(lfc_meth < 0) %>%
  mutate(pos=paste0(chr, ":", start, "-", end)) %>%
  mutate(category_2 = ifelse(meth_g2 <= 0.5*meth_g1, "yes", "no")) %>%
  group_by(cmp, category_2) %>%
  reframe(n=n())

