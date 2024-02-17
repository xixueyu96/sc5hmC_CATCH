##library packages
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggridges)

## load input
chr_windows <- readRDS("data/processed/mm10.CpG_density.win1k.noN.rds")

df_5hmC <-
  fread(
    "data/processed/mouse_embryo.5hmC.1X50.aim.1kb.230104.bed",
    header = T,
    stringsAsFactors = F
  )
df_5hmC_5mC_nor <-
  fread(
    "data/processed/mouse_embryo.5hmC_5mC_nor.1X50.aim.1kb.230104.bed",
    header = T,
    stringsAsFactors = F
  )

df_5hmC %<>% filter(end %in% df_5hmC_5mC_nor$end)
df_5hmC_5mC_nor %<>% mutate(start = start + 1)

stage_list <- colnames(df_5hmC) %>% .[4:length(.)]
colnames(df_5hmC_5mC_nor) <- c("chrom", "start", "end", stage_list)

## plot 5hmC

df_5hmC %<>%
  inner_join(chr_windows %>% as.data.table(),
             by=c("chrom"="seqnames",
                  "start"="start","end"="end")) %>%
  mutate(
    cpg_density = case_when(
      cpg_density >= 0.08 ~ 0.08,
      cpg_density <= 0.03 ~ 0.03,
      TRUE ~ cpg_density)) %>%
  mutate(cpg_density_bin = cut(
    cpg_density,
    breaks = seq(0.03, 0.08, 0.01),
    include.lowest = T))
# hist(df_5hmC$`2C` %>% .[.<0.1], breaks = 100)

df_5hmC_long <- df_5hmC %>%
  dplyr::select(.i_window, cpg_density_bin, stage_list) %>%
  reshape2::melt(id.var=c(".i_window", "cpg_density_bin")) %>%
  filter(!is.na(value))

df_5hmC_long$variable <- factor(
  df_5hmC_long$variable,
  levels = c(
    "Oocyte", "Sperm", "Zygote",
    "2C", "4C", "8C",
    "Morula", "Blast")
)

df_5hmC_long %>%
  dplyr::filter(value <=as.numeric(quantile(df_5hmC_long$value, probs = c(0.99)))) %>%
  dplyr::filter(value > 0) %>%
  ggplot(aes(x = value, y = cpg_density_bin, group=cpg_density_bin)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01,fill="#80B1D3") +
  facet_wrap(.~variable, scales = "free_y") +
  labs(y="CpG density",x="5hmCpG level") +
  theme_bw(base_size = 15)+
  theme_ridges()+
  theme(strip.background = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(vjust = 0))

ggsave(
  "viz/Fig3C.5hmC.240217.pdf",
  width = 27.5,
  height = 24,
  units = "cm"
)

## plot 5hmC_5mC_nor

df_5hmC_5mC_nor %<>%
  inner_join(
    chr_windows %>% as.data.table,
    by=c("chrom"="seqnames",
         "start"="start","end"="end")) %>%
  mutate(
    cpg_density = case_when(
      cpg_density >= 0.08 ~ 0.08,
      cpg_density <= 0.03 ~ 0.03,
      TRUE ~ cpg_density)) %>%
  mutate(cpg_density_bin = cut(
    cpg_density,
    breaks = seq(0, 0.08, 0.01),
    include.lowest = T))

df_5hmC_5mC_nor_long <- df_5hmC_5mC_nor %>%
  dplyr::select(.i_window, cpg_density_bin, stage_list) %>%
  reshape2::melt(id.var=c(".i_window", "cpg_density_bin")) %>%
  filter(!is.na(value))

df_5hmC_5mC_nor_long$variable <- factor(
  df_5hmC_5mC_nor_long$variable,
  levels = c(
    "Oocyte", "Sperm", "Zygote",
    "2C", "4C", "8C",
    "Morula", "Blast")
)

df_5hmC_5mC_nor_long %>%
  dplyr::filter(value <=as.numeric(quantile(df_5hmC_5mC_nor_long$value, probs = c(0.99)))) %>%
  dplyr::filter(value > 0) %>%
  # ggplot(aes(x=value, y=cpg_density)) +
  ggplot(aes(x = value/100, y = cpg_density_bin, group=cpg_density_bin)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01,fill="#80B1D3") +
  facet_wrap(.~variable, scales = "free_y") +
  labs(y="CpG density",x="5hmCpG level (normed by 5mCpG)") +
  theme_bw(base_size = 15)+
  theme_ridges()+
  theme(strip.background = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(vjust = 0))

ggsave(
  "viz/Fig3C.5hmC_5mC_nor.240217.pdf",
  width = 27.5,
  height = 24,
  units = "cm"
)

## stat
cor_test_with_CpG_density <- function(x, dataset){

  upper_limit <- quantile(
    dataset %>% select(x) %>%
      filter(!is.na(get(x))) %>% pull,
    probs = 0.99)

  tmp <- dataset %>%
    filter(get(x) <= upper_limit & get(x) > 0) %>%
    # filter(get(x) <= upper_limit) %>%
    dplyr::select(x, cpg_density) %>%
    as.data.frame()

  cor_res <- cor.test(tmp[, 1], tmp[, 2], method="spearman")

  return(
    data.frame(stage=x, cor=cor_res$estimate, pvalue=cor_res$p.value)
  )
}

cor_res_5hmC <- do.call(
  "rbind",
  lapply(
    stage_list,
    cor_test_with_CpG_density,
    df_5hmC)
)

cor_res_5hmC

write.table(
  cor_res_5hmC,
  file = "publish/Fig3C.5hmC_cor_test.240217.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F
)

cor_res_5hmC_5mC_nor <- do.call(
  "rbind",
  lapply(
    stage_list,
    cor_test_with_CpG_density,
    df_5hmC_5mC_nor)
)

cor_res_5hmC_5mC_nor

write.table(
  cor_res_5hmC,
  file = "publish/Fig3C.5hmC_5mC_nor.cor_test.240217.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F
)

