library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
# library(ggpubr)

gene_list <- c("Tet1", "Tet2", "Tet3", "Tdg", "Dnmt1", "Uhrf1", "Dnmt3a", "Dnmt3b", "Dnmt3l")

if(F) {
  mouse_preim <- readRDS("data/processed/Mouse_preimplan.seurat.rds")

  stage_level <- c("M2","ZY","E2C", "L2C", "4C", "8C", "Morula", "ICM","TE")

  mouse_preim$Stage <- factor(mouse_preim$Stage, levels = stage_level)

  VlnPlot(mouse_preim, features = gene_list, group.by = "Stage")

  ## FigS3A by SCAN-seq
  df3 <- mouse_preim@assays$RNA@data %>%
    as.matrix() %>% .[gene_list, ] %>%
    reshape2::melt() %>%
    left_join(
      mouse_preim@meta.data[,c("Cell_id", "Stage")],
      by=c("Var2"="Cell_id")
    )
  df3$Stage <- factor(df3$Stage, levels = stage_level)

  todo_index <- 1:length(gene_list)
  barplot_fun <- function(x){
    g <- gene_list[x]
    p <- df3 %>% filter(Var1==g) %>%
      ggbarplot(
        x = "Stage", y = "value",fill="#6967aa",
        add = "mean_se", error.plot = "upper_errorbar")+
      ylim(min(df3$value), max(df3$value))+
      labs(y="Expression level", title=g)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.title = element_text(hjust = 0.5, face = "italic"))
    return(p)
  }

  p_list <- lapply(todo_index, barplot_fun)
  do.call("plot_grid", p_list)

  ggsave(
    "viz/FigS3A.5mC_related_gene_expression.pdf",
    width = 7.5, height = 5.7
  )
}

## FigS3A by Qiao Jie

if(F) {
  rna_mtx <- readRDS('data/processed/GSE136714.aggregated_gene_matrix.Rds')

  ### normalization (count -> FPKM) ###
  gene_anno <- readRDS('data/public/gencode.vM11.gene_length.Rds')

  ## check gene symbol & ensembl
  gene_anno %>%
    filter(symbol %in% gene_list) %>%
    reframe(duplicated(symbol))

  rna_mtx <- t(t(rna_mtx)/colSums(rna_mtx)*1e6)

  gene_length <- gene_anno %>%
    filter(symbol %in% gene_list) %>%
    dplyr::select(exonic.gene.sizes, symbol)
  rownames(gene_length) <- gene_length$symbol

  plot_mtx <- rna_mtx[gene_list,]/gene_length[gene_list,"exonic.gene.sizes"]*1e3

  ### generate anno for cell ###
  plot_df <- plot_mtx %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    reshape2::melt("gene") %>%
    filter(!grepl("Mixed", variable)) %>%
    mutate(variable = gsub("_embryo|Late", "", variable)) %>%
    tidyr::separate(
      col = "variable",
      into = c("stage_embryo", "cell_id"),
      sep = "_") %>%
    tidyr::separate(
      col = "stage_embryo",
      into = c("stage", "embryo"),
      sep = "(?<=\\D)(?=\\d)")

  plot_df$stage <- factor(plot_df$stage, levels = c("Zygote", "2cell", "4cell", "8cell", "16cell", "32cell"))

  todo_index <- 1:length(gene_list)
  barplot_fun <- function(x){
    g <- gene_list[x]

    if(grepl("Tet", g)){
      p <- plot_df %>% filter(gene==g) %>%
        ggbarplot(
          x = "stage", y = "value",fill="#6967aa",
          add = "mean_se", error.plot = "upper_errorbar")+
        ylim(min(plot_df$value), max(plot_df$value))+
        labs(y="Expression level", title=g)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              plot.title = element_text(hjust = 0.5, face = "italic"))
      return(p)
    } else {
      p <- plot_df %>% filter(gene==g) %>%
        ggbarplot(
          x = "stage", y = "value",fill="#6967aa",
          add = "mean_se", error.plot = "upper_errorbar")+
        # ylim(min(plot_df$value), max(plot_df$value))+
        labs(y="Expression level", title=g)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              plot.title = element_text(hjust = 0.5, face = "italic"))
      return(p)
    }
  }

  p_list <- lapply(todo_index, barplot_fun)
  do.call("plot_grid", p_list)

  ggsave(
    "viz/FigS3A.5mC_related_gene_expression.240509.QJ.v2.pdf",
    width = 7.5, height = 5.7
  )

  plot_df %>%
    group_by(gene, stage) %>%
    filter(grepl("Tet", gene)) %>%
    summarise(mean=mean(value))
}


