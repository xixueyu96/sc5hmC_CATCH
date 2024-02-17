## library packages
library(dplyr)
library(ggplot2)
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

## plot result
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
