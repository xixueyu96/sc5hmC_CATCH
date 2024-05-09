library(dplyr)
library(data.table)
library(GenomicRanges)

total_5hmC <- fread("data/source/TET123_5hmC_site.bed") %>%
  setNames(c("seqnames", "start", "end", "base", "cov", "met")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  .[seqnames(.) %in% paste0("chr", c(1:19, "X", "Y", "MT")),]

TET3KO_5hmC <- fread("data/source/TET3cKO_EZ.5hmC_site.bed") %>%
  setNames(c("seqnames", "start", "end", "met",  "cov")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  .[seqnames(.) %in% paste0("chr", c(1:19, "X", "Y", "MT")),]

TET3KO_non_5hmC <- fread("data/source/TET3cKO_non_5hmC_site.bed") %>%
  setNames(c("seqnames", "start", "end", "met",  "cov")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  .[seqnames(.) %in% paste0("chr", c(1:19, "X", "Y", "MT")),]

oocyte_5hmC <- fread("data/source/Oocyte.5hmC_site.bed") %>%
  setNames(c("seqnames", "start", "end", "base", "cov",  "meth")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  .[seqnames(.) %in% paste0("chr", c(1:19, "X", "Y", "MT")),]

library(ChIPpeakAnno)
res <- makeVennDiagram(
  Peaks=list(TET3KO_5hmC, oocyte_5hmC),
  NameOfPeaks=c("TET3cKO_5hmC", "Oocyte_5hmC")
)

res <- makeVennDiagram(
  Peaks=list(oocyte_5hmC, total_5hmC),
  NameOfPeaks=c("Oocyte_5hmC", "total_5hmC")
)

# devtools::install_github("js229/Vennerable")
library(Vennerable)
venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}
v <- venn_cnt2venn(res$vennCounts)
plot(v)

TET12_5hmC <- overlapsAny(total_5hmC, TET3KO_5hmC)
TET3_5hmC <- overlapsAny(total_5hmC, TET3KO_non_5hmC)

TET12_5hmC <- total_5hmC[TET12_5hmC]
TET3_5hmC <- total_5hmC[TET3_5hmC]

saveRDS(TET12_5hmC, file = "data/processed/TET12_specific_5hmC_site.rds")
saveRDS(TET3_5hmC, file = "data/processed/TET3_specific_5hmC_site.rds")

