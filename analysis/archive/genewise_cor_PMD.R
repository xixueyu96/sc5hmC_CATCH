library(data.table)

gene_cor_res <- read.table(
  "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
  header = T, stringsAsFactors = F
)

pos_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho > 0) %>%
  arrange(rho)

neg_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho < 0) %>%
  arrange(-rho)

gene_gr <- fread(
  "data/public/mm10.genebody.bed",
  header = T, stringsAsFactors = F) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

pos_genes_gr <- gene_gr[gene_gr$gene %in% pos_genes$gene]
neg_genes_gr <- gene_gr[gene_gr$gene %in% neg_genes$gene]

## download from https://zwdzwd.github.io/pmd
pmd_gr <- fread(
  "data/public/PMD_coordinates_mm10.bed.gz",
  header = F, stringsAsFactors = F) %>%
  setNames(c("chr", "start", "end", "score", "type", "iscommonPMD")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

pmd_gr <- pmd_gr[!is.na(pmd_gr$type) & pmd_gr$type == "PMD"]

pos_hit <- findOverlaps(pos_genes_gr, pmd_gr)
pos_hit <- pos_genes_gr[pos_hit@from]
length(pos_hit)/length(pos_genes_gr)

neg_hit <- findOverlaps(neg_genes_gr, pmd_gr)
neg_hit <- neg_genes_gr[neg_hit@from]
length(neg_hit)/length(neg_genes_gr)

fisher.test(
  data.frame(
    "pos" = c(length(pos_hit), length(pos_genes_gr)),
    "neg" = c(length(neg_hit), length(neg_genes_gr))
  )
)

overlaps <- intersect(neg_genes_gr, pmd_gr)
total_overlap_bases <- sum(width(overlaps))
total_overlap_bases/sum(width(neg_genes_gr))


