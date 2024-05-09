## library packages
library(dplyr)
library(magrittr)

## load input
mc_mtx <- read.table(
  "data/processed/mC.genebody_mtx.tsv",
  header = T, stringsAsFactors = F
)
hmc_mtx <- read.table(
  "data/processed/hmC.genebody_mtx.tsv",
  header = T, stringsAsFactors = F
)
rna_mtx <- read.table(
  "data/processed/RNA.genebody_mtx.tsv",
  header = T, stringsAsFactors = F
)

stage_to_test <-
  c("Two.cell",
    "Four.cell",
    "Eight.cell",
    "Blast",
    "Morula",
    "Zygote")

hmc_mtx %<>% tibble::column_to_rownames("gene") %>%
  dplyr::select(all_of(stage_to_test))
mc_mtx %<>% tibble::column_to_rownames("gene") %>%
  dplyr::select(all_of(stage_to_test))
rna_mtx %<>% tibble::column_to_rownames("gene") %>%
  dplyr::select(all_of(stage_to_test))

## calc correlation
gene_cor <- function(i){
  # print(i)
  result <- try(
    res <- cor.test(
      hmc_mtx %>% .[i,] %>% unlist,
      mc_mtx %>% .[i,] %>% unlist,
      method = "spearman"),
    silent = T
  )
  if(class(result)=="try-error"){
    return(NULL)
  }else {
    df <- data.frame(
      gene = i,
      rho = res$estimate,
      pval = res$p.value
    )
    return(df)
  }
  gc()
}

gene_cor_res_list <- lapply(
  rownames(hmc_mtx), gene_cor
)

gene_cor_res <- do.call(
  rbind, gene_cor_res_list
)

gene_cor_res %<>%
  # mutate(padj = p.adjust(pval, method = "fdr")) %>%
  mutate(
    group = case_when(
      pval < 0.05 & rho > 0 ~"pos",
      pval < 0.05 & rho < 0 ~ "neg",
      TRUE ~ "none"
      )
  )
table(gene_cor_res$group)

write.table(
  gene_cor_res,
  file = "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F
)

## metadata of genes
library(GenomicRanges)

gene_df <- read.table(
  "data/public/mm10.genebody.bed",
  header = T, stringsAsFactors = F
)

gene_gr <- makeGRangesFromDataFrame(gene_df, keep.extra.columns = T)
gene_gr <- sort(gene_gr)
gene_gr

## all the CpGs
all_cpg <- readRDS("data/handmade/mm10.CpG.rds")
all_cpg <- sort(all_cpg)

## from jzh @biubiu-1s
hits <- findOverlaps(all_cpg, gene_gr)
hits_rle <- rle(sort(hits@to))

gene_gr$cpg_density <- 0
gene_gr$cpg_density[hits_rle$values] <- hits_rle$lengths/width(gene_gr)[hits_rle$values]

saveRDS(gene_gr, file = "data/handmade/mm10.genebody.rds")

gc()

