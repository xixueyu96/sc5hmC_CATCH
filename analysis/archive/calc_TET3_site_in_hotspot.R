library(GenomicRanges)
library(data.table)
library(ChIPpeakAnno)
library(ggpubr)

region_pal <- c("hotspot"="#7ca9d5", "matched"="#e7df68", "other"="#e5545b")

TET3_5hmC <- readRDS("data/processed/TET3_specific_5hmC_site.rds")
# TET12_5hmC <- readRDS("data/processed/TET12_specific_5hmC_site.rds")
merged_gr <- readRDS("data/processed/hotspot.5hmC.1X50.240218.rds")
total_5hmC <- fread("data/source/TET123_5hmC_site.bed") %>%
  setNames(c("seqnames", "start", "end", "base", "cov", "met")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  .[seqnames(.) %in% paste0("chr", c(1:19, "X", "Y", "MT")),]


## get overlapping count

get_cnt <- function(gr){

  cnt <- lapply(
    c("other", "matched", "hotspot"),

    function(x){
      hit <- findOverlaps(
        merged_gr %>% subset(region==x), gr)
      return(
        data.frame(
          # "region"=x,
          "count"=length(hit@to),
          row.names = x
      ))}
  )

  cnt_df <- do.call(rbind, cnt) %>%
    tibble::rownames_to_column("region") %>%
    mutate(per = count/sum(count)*100)

  return(cnt_df)
}

## TET3
tet3_cnt_df <- get_cnt(TET3_5hmC)

labs <- paste0(round(tet3_cnt_df$per), "%")

ggdonutchart(
  tet3_cnt_df, "count",
   label = labs,
   lab.pos = "in", lab.font = "black",
   fill = "region", color = "white",
   palette = region_pal
)

ggsave(
  "viz/TET3_specific_5hmC_in_hotspot.dounut.pdf",
  width = 5, height = 5
)

## total
total_cnt_df <- get_cnt(total_5hmC)

labs <- paste0(round(total_cnt_df$per), "%")

ggdonutchart(
  total_cnt_df, "count",
  label = labs,
  lab.pos = "in", lab.font = "black",
  fill = "region", color = "white",
  palette = region_pal
)

ggsave(
  "viz/total_5hmC_in_hotspot.dounut.pdf",
  width = 5, height = 5
)

## test
dat <- data.frame(
  "total" = total_cnt_df$count,
  "tet3" = tet3_cnt_df$count,
  row.names = total_cnt_df$region) %>% t

pdf("viz/TET3_specific_5hmC_in_hotspot.pdf", width = 4, height = 8)
mosaicplot(
  dat[,3:1],
  main = "Mosaic plot",
  color = TRUE
)
dev.off()

fisher.test(dat[2:1,3:2])

# p-value = 0.02797
# odds ratio = 2.259869

fisher.test(dat[2:1,2:1])

fisher.test(dat[2:1,c(3,1)])

