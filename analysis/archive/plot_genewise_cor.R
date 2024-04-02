library(dplyr)
library(ggplot2)

gene_cor_res <- read.table(
  "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
  header = T, stringsAsFactors = F
)

## volcano
gene_cor_res %>%
  ggplot(aes(x=rho, y=-log10(pval)))+
  geom_point(aes(color=group), shape=21, size=3) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash") +
  geom_vline(xintercept = 0, linetype="longdash") +
  scale_color_manual(values = c("pos"="#e51d17", "neg"="#3d51a2", "none"="grey30")) +
  theme_presentation() +
  theme(aspect.ratio = 1)

ggsave(
  "viz/Gene_wise_cor.5hmC_5mC.volcano.240217.pdf",
  width = 6, height = 6.5
)

## pos & neg vs 5hmC, 5mC, RNA
stage_pal <- c(
  "Zygote" = "#5092c6", "Two.cell"="#f47f70", "Four.cell"="#f7b463",
  "Eight.cell" = "#b5df6e", "Morula"="#f9cce4", "Blast"="#d9d9d9",
  "Oocyte" = "#94d3c8", "Sperm" = "#da9dfd"
)

pos_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho > 0) %>%
  pull(gene)
neg_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho < 0) %>%
  pull(gene)

plot_MO_genes <- function(genes){

  plot_df <- rbind(
    hmc_mtx %>% .[genes,] %>% reshape2::melt() %>%
      mutate(omic="5hmC") %>% filter(value < 0.03),
    rna_mtx %>% .[genes,] %>% reshape2::melt() %>%
      mutate(omic="RNA") %>% filter(value < 0.03),
    mc_mtx %>% .[genes,] %>% reshape2::melt() %>%
      mutate(omic="5mC") %>% filter(value < 0.3)
  )
  plot_df$variable <- factor(
    plot_df$variable,
    levels =
      c("gene",
        "Oocyte",
        "Sperm",
        "Zygote",
        "Two.cell",
        "Four.cell",
        "Eight.cell",
        "Morula",
        "Blast")
  )

  p1 <- plot_df %>%
    # filter(variable %in% colnames(hmc_mtx)[3:ncol(hmc_mtx)]) %>%
    ggplot(aes(x=variable, y=value))+
    # ggbeeswarm::geom_quasirandom(aes(color=variable)) +
    geom_violin(aes(color=variable), scale = "width", size=1) +
    geom_boxplot(aes(color=variable), outlier.shape = NA, width=.3, size=1) +
    facet_wrap(.~omic, scales = "free", nrow=1) +
    scale_color_manual(values = stage_pal) +
    labs(y="Modification/Expression") +
    theme_presentation() +
    theme(aspect.ratio = 1,
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())

  # p1
  return(p1)
}

p1 <- plot_MO_genes(pos_genes)
p2 <- plot_MO_genes(neg_genes)

cowplot::plot_grid(p1, p2, ncol = 1)

ggsave(
  "viz/Gene_wise_cor.5hmC_5mC.20240217.pdf",
  width = 12, height = 10
)

## gene length
gene_gr <- readRDS("data/handmade/mm10.genebody.rds")

plot_df <- inner_join(
  gene_cor_res,
  gene_gr %>% as.data.frame(),
  by = c("gene")
)

plot_df %>%
  dplyr::select(group, cpg_density, width) %>%
  reshape2::melt(id.var="group") %>%
  filter(group!="none") %>%
  ggplot(aes(x=group, y=log10(value)))+
  geom_violin(aes(color=group)) +
  geom_boxplot(aes(color=group), outlier.shape = NA, size=1, width=.3) +
  scale_color_manual(values = c("pos"="#e51d17", "neg"="#3d51a2", "none"="grey30")) +
  ggpubr::stat_compare_means(comparisons = list(c("pos", "neg"))) +
  facet_wrap(.~variable, scale="free_y") +
  egg::theme_presentation()+
  theme(aspect.ratio = 1.2)

ggsave(
  "viz/Gene_wise_cor.CpGdensity_length.240217.pdf",
  width = 10, height = 5
)

## H3K36me3
k36_mtx <- lapply(
  list.files("data/source/h3k36me3_genebody/", full.names = T),
  function(x){
    df <- read.table(x, header = F, stringsAsFactors = F)
    sp <- strsplit(basename(x), "[.]")[[1]][1]
    colnames(df) <- c("chr", "start", "end", "gene", sp)
    df <- df %>% .[,c("gene", sp)]
    print(sapply(df, class))
    return(df)
  }
)

k36_mtx <-
  k36_mtx %>% purrr::reduce(function(x, y)
    base::merge(x, y, by = c("gene"))) %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.))

plot_df <- inner_join(
  k36_mtx, gene_cor_res, by="gene") %>%
  dplyr::select(group, `1C`, `Late_2C`, `8C`, `ICM`) %>%
  reshape2::melt(id.var="group") %>%
  filter(group!="none")

stat.test <- plot_df %>%
  group_by(variable) %>%
  summarise(pval=wilcox.test(value~group)$p.value)
stat.test

plot_df %>%
  ggplot(aes(x=variable, y=value))+
  geom_boxplot(aes(color=group), outlier.shape = 21, size=1) +
  # geom_hline(yintercept = 0, linetype="longdash") +
  scale_color_manual(values = c("pos"="#e51d17", "neg"="#3d51a2", "none"="grey30")) +
  # stat_compare_means(comparisons = list(c("pos", "neg"))) +
  # stat_compare_means(aes(group = variable)) +
  # ggsignif::geom_signif(comparisons = list(c("pos", "neg"))) +
  annotate(
    "text", label = format(stat.test$pval, digit=2),
    x = stat.test$variable, y = 4.8) +
  labs(y="Averaged signals", title = "H3K36me3")+
  egg::theme_presentation()+
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = .5),
        axis.title.x = element_blank())

ggsave(
  "viz/Gene_wise_cor.H3K36me3.violin.240217.pdf",
  width = 8, height = 8
)


## H3K27me3
k27_mtx <- lapply(
  list.files("data/source/H3K27me3_genebody/", full.names = T),
  function(x){
    df <- read.table(x, header = F, stringsAsFactors = F)
    sp <- strsplit(basename(x), "[.]")[[1]][1]
    # sp <- strsplit(sp, "_")[[1]][1]
    colnames(df) <- c("chr", "start", "end", "gene", sp)
    df <- df %>% .[,c("gene", sp)]
    print(sapply(df, class))
    return(df)
  }
)

k27_mtx <-
  k27_mtx %>% purrr::reduce(function(x, y)
    base::merge(x, y, by = c("gene"))) %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.)) %>%
  reshape2::melt(id.var="gene") %>%
  mutate(variable=stringr::str_split(variable, "_",simplify = T)[,2]) %>%
  group_by(gene, variable) %>%
  summarise(value_avg=mean(value)) %>%
  ungroup() %>%
  tidyr::spread(value = value_avg, key = variable)

plot_df <- inner_join(
  k27_mtx, gene_cor_res, by="gene") %>%
  dplyr::select(-rho, -pval, -gene) %>%
  reshape2::melt(id.var="group") %>%
  filter(group!="none")

stage_ordre <- c("2cell", "4cell", "8cell", 'morula', 'ICM', 'TE')

stat.test <- plot_df %>%
  group_by(variable) %>%
  summarise(pval=wilcox.test(value~group)$p.value) %>%
  arrange(variable, levels = stage_ordre)
stat.test

plot_df$variable <- factor(plot_df$variable, levels = stage_ordre)

plot_df %>%
  ggplot(aes(x=variable, y=value))+
  geom_boxplot(aes(color=group), outlier.shape = 21, size=1) +
  # geom_hline(yintercept = 0, linetype="longdash") +
  scale_color_manual(values = c("pos"="#e51d17", "neg"="#3d51a2", "none"="grey30")) +
  # stat_compare_means(comparisons = list(c("pos", "neg"))) +
  # stat_compare_means(aes(group = variable)) +
  # ggsignif::geom_signif(comparisons = list(c("pos", "neg"))) +
  ylim(c(0,0.4)) +
  annotate(
    "text", label = format(stat.test$pval, digit=2),
    x = stat.test$variable, y = 0.37) +
  labs(y="Averaged signals", title = "H3K27me3")+
  egg::theme_presentation()+
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = .5),
        axis.title.x = element_blank())

ggsave(
  "viz/Gene_wise_cor.H3K27me3.violin.240217.pdf",
  width = 8, height = 8
)

## cross-omics correlation
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

plot_df <- inner_join(
  mc_mtx, hmc_mtx,
  by = c("gene"),
  suffix = c("_mC", "_hmC")) %>%
  inner_join(
    rna_mtx %>%
      rename_with(
        .fn = function(.x) {
          paste0(.x, "_rna")},
        .cols = -c(gene)),
    by = c("gene" = "gene")
  )

cor_matrix <- plot_df %>%
  filter(gene %in% neg_genes) %>%
  dplyr::select(starts_with("Two")) %>%
  filter(`Two.cell_rna` > 0) %>%
  cor(method = "spearman")

corrplot::corrplot(
  cor_matrix,
  type = "lower",
  method = "square",
  addCoef.col = "black",
  tl.col = "black", tl.srt = 45
)

## enrichedheatmap
library(EnrichedHeatmap)
library(rtracklayer)
library(data.table)

genes <- fread("data/public/mm10.gencode.p4.protein_coding.bed") %>%
  setNames(c("seqnames", "start", "end", "ensembl", "symbol", "strand")) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  filter(symbol %in% c(pos_genes, neg_genes)) %>%
  mutate(
    cor = case_when(
      symbol %in% pos_genes ~ "pos",
      symbol %in% neg_genes ~ "neg"),
    seqnames = paste0("chr", seqnames)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  .[seqnames(.) %in% paste0("chr", c(1:19, "X", "Y", "MT")),]
genes

tss <- promoters(genes, upstream = 0, downstream = 1)
tss[1:5]

line_col <- c("pos"="#e51d17", "neg"="#3d51a2")

rna_mtx <- fread("data/public/GSE71434_FPKM_stage.txt.gz", header = T) %>%
  filter(`#Gene` %in% c(pos_genes, neg_genes)) %>%
  mutate(
    cor = case_when(
      `#Gene` %in% pos_genes ~ "pos",
      `#Gene` %in% neg_genes ~ "neg"))

Heatmap(
  rna_mtx %>% select(PN5_zygote, `2cell_early`,`2cell_late`),
  # col = c("white", "darkgreen"),
  name = "PolII_E",
  # top_annotation = top_anno,
  row_split = rna_mtx$cor,
  column_title = "PolII in PN5",
  width=unit(4,"cm"),
  cluster_rows = F,
  cluster_columns = F,
  # use_raster=T,
  # raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)

top_anno <- HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = rev(line_col))))


top_anno_pol <- HeatmapAnnotation(
  lines = anno_enriched(
    gp = gpar(col = rev(line_col)),
    ylim = c(0,8),
    axis_param = list(
    side = "right",
    at = seq(2,8,2),
    labels = seq(2,8,2))
  )
)

ref <- import("data/public/GSM1845268_2cell_late_K4me3_rep1.bed")
ref$score <- as.numeric(ref$name)

mat0 <-  normalizeToMatrix(
  ref, tss, value_column = "score",
  extend = 3000, mean_mode = "absolute", w = 50,
  keep = c(0, 0.99)
)
mat0

# rm(ref)

ht0 <- EnrichedHeatmap(
  mat0, col = c("white", "orange"),
  name = "H3K4me3",
  top_annotation = top_anno,
  row_split = tss$cor,
  column_title = "H3K4me3(L2C)",
  width=unit(4,"cm"),
  use_raster=T,
  raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)
ht0

pol2_E2C <- import("data/public/GSM4010568_PN5_RP_rep1.mm10.bw")
pol2_E2C

mat1 <-  normalizeToMatrix(
  pol2_E2C, tss, value_column = "score",
  extend = 3000, mean_mode = "w0", w = 50,
  keep = c(0, 0.99)
)
mat1

pol2_E_col_fun <- circlize::colorRamp2(
  c(0, max(mat1)), c("white", "darkgreen")
)

rm(pol2_E2C)

ht1 <- EnrichedHeatmap(
  mat1, col = pol2_E_col_fun,
  name = "PolII_ZY",
  top_annotation = top_anno_pol,
  row_split = tss$cor,
  column_title = "PolII(PN5)",
  width=unit(4,"cm"),
  use_raster=T,
  raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)
ht1

pol2_L2C <- import("data/public/GSM4010570_2cell_early_RP_rep1.mm10.bw")
pol2_L2C

mat2 <-  normalizeToMatrix(
  pol2_L2C, tss, value_column = "score",
  extend = 3000, mean_mode = "w0", w = 50,
  keep = c(0, 0.99)
)
mat2
rm(pol2_L2C)

ht2 <- EnrichedHeatmap(
  mat2, col = pol2_E_col_fun,
  name = "PolII_L2C",
  top_annotation = top_anno_pol,
  row_split = tss$cor,
  column_title = "PolII(L2C)",
  width=unit(4,"cm"),
  use_raster=T,
  raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)
ht1 + ht2

h3k36me3 <- import("data/public/Late_2C_H3K36me3.bw")
h3k36me3

mat4 <-  normalizeToMatrix(
  h3k36me3, genes, value_column = "score",
  extend = 3000, mean_mode = "w0", w = 50,
  keep = c(0, 0.99)
)
mat4
rm(h3k36me3)

ht4 <- EnrichedHeatmap(
  mat4, col = c("white", "blue"),
  name = "H3K36me3",
  top_annotation = top_anno,
  row_split = tss$cor,
  column_title = "H3K36me3(L2C)",
  width=unit(4,"cm"),
  use_raster=T,
  raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)

h3k27me3 <- import("data/public/GSM2082673_2cell.H3K27me3.1.bedGraph")
h3k27me3

mat5 <-  normalizeToMatrix(
  h3k27me3, genes, value_column = "score",
  extend = 3000, mean_mode = "w0", w = 50,
  keep = c(0, 0.99)
)
mat5
rm(h3k27me3)

ht5 <- EnrichedHeatmap(
  mat5, col = c("white", "red"),
  name = "H3K27me3",
  top_annotation = top_anno,
  row_split = tss$cor,
  column_title = "H3K27me3(2C)",
  width=unit(4,"cm"),
  use_raster=T,
  raster_resize_mat=mean,
  raster_quality = 3,
  heatmap_legend_param = list(
    # direction = "horizontal",
    border = "black",
    legend_height = unit(2, "cm")
  )
)

lgd <- Legend(at = names(line_col), title = "Genes",
             type = "lines", legend_gp = gpar(col = line_col))

# cgi <- import("data/public/mm10.CGI.bed.gz")
#
# mat_cgi <- normalizeToMatrix(cgi, tss, mean_mode = "absolute")
#
# ht <- EnrichedHeatmap(
#   mat_cgi, col = c("white", "orange"), name = "CGI",
#   top_annotation = top_anno,
#   row_split = tss$cor,
#   column_title = "H3K4me3 in 2C",
#   width=unit(5,"cm")
# )

pdf("viz/Gene_wise_cor.heatmap.240402.pdf",
    width = 12, height = 6)
draw(
  ht1 + ht2 + ht0 + ht4 + ht5,
  annotation_legend_list = list(lgd),
  ht_gap = unit(c(10,10,10,10), "mm"),
  heatmap_legend_side = "right"
)
dev.off()
