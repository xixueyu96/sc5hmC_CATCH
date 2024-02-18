library(dplyr)
library(GenomicRanges)
library(MatchIt)
library(data.table)

## load input
hotspot_gr <- LOLA::readBed("data/source/1111_hmC_ratio_1X50.merge.1k.v5s4.hotspot.bed")
seqlevels(hotspot_gr) <- hotspot_gr %>%
  seqnames %>% unique %>% stringr::str_sort(numeric = T)
hotspot_gr <- sort(hotspot_gr)

cpg_density <- readRDS("data/handmade/mm10.CpG_density.win1k.noN.rds")
cpg_density

seqlevels(cpg_density)
cpg_density <- sort(cpg_density)

## add CpG density
hits <- findOverlaps(hotspot_gr, cpg_density)
cpg_density$hotspot <- 0
cpg_density$hotspot[hits@to] <- 1
table(cpg_density$hotspot)

## sample non-hotspot regions
non_hotspot_gr <- cpg_density %>% .[.$hotspot==0]
non_hotspot_gr <- non_hotspot_gr[sample(length(non_hotspot_gr), size = 5e4)]

sampled_gr <- c(
  cpg_density %>% .[.$hotspot==1],
  non_hotspot_gr
)
seqlevels(sampled_gr)
sampled_gr <- sort(sampled_gr)

gc()

## stat hotspot against cpg_density
wilcox.test(
  sampled_gr %>% .[.$hotspot == 1] %>% .$cpg_density,
  sampled_gr %>% .[.$hotspot == 0] %>% .$cpg_density,
  alternative = "less"
) ## p-value = 6.693e-06

## matching
m.out1 <- matchit(
  hotspot ~ cpg_density,
  data = sampled_gr %>% as.data.table(),
  method = "nearest",
  distance = "glm", link = "probit"
)
summary(m.out1)

gc()

plot(
  m.out1, type = "jitter",
  interactive = FALSE
)
plot(
  m.out1, type = "density",
  interactive = FALSE,
  which.xs = ~ cpg_density
)

wilcox.test(
  get_matches(m.out1) %>% subset(hotspot==1) %>% .$cpg_density,
  get_matches(m.out1) %>% subset(hotspot==0) %>% .$cpg_density,
  alternative = "less"
) # p-value = 0.5

control_gr <- get_matches(m.out1) %>%
  subset(hotspot==0) %>%
  makeGRangesFromDataFrame()
control_gr <- sort(control_gr)

saveRDS(control_gr, file = "data/processed/hotspot.cpg_density_matched.rds")

