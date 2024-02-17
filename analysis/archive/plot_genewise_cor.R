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



