library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

## Fig3E by SCAN-seq
if(F) {
  gse_path <- "/mnt/f/Project/sc5hmC/data/ref/total_RNA/expr_mtx/mm10/repeat_GSR/"

  merge2 <- function(x,y)base::merge(x,y,by=c("Geneid"),all=T)

  rna_mtx <- lapply(
    list.files(gse_path, full.names = T),
    function(x){
      df <- read.table(x, header = T) %>% dplyr::select(-Length)
      colnames(df)[2] <- strsplit(colnames(df)[2], "[.]")[[1]][1]
      df[,2] <- df[,2]/sum(df[,2])*1e4
      return(df)}) %>%
    purrr::reduce(merge2)

  ## names repeats in family L1, L2, Alu and MIR
  repeat_name <- readRDS("data/processed/mm10.repeat_name.rds")

  plot_mtx <- rna_mtx %>%
    mutate(
      repeat_class = case_when(
        Geneid %in% repeat_name$L1 ~ "L1",
        Geneid %in% repeat_name$L2 ~ "L2",
        Geneid %in% repeat_name$Alu ~ "Alu",
        Geneid %in% repeat_name$MIR ~ "MIR",
        Geneid %in% repeat_name$ERVL ~ "ERVL",
        Geneid %in% repeat_name$ERVK ~ "ERVK",
        TRUE ~ "Others")) %>%
    filter(repeat_class!="Others") %>%
    reshape2::melt(
      id.vars = c("Geneid", "repeat_class"),
      variable.name = "sample",
      value.name = "count") %>%
    tidyr::separate(
      col = "sample",
      into = c("srr", "stage"),
      sep = "_")

  stage_level <- c("2C", "4C", "8C", "Morula", "ICM", "TE")
  repeat_level <- c('L1', 'L2', 'Alu', 'MIR', 'ERVL', 'ERVK')

  plot_mtx$stage <- factor(plot_mtx$stage, levels = stage_level)
  plot_mtx$repeat_class <- factor(plot_mtx$repeat_class, levels = repeat_level)

  plot_mtx %>%
    group_by(repeat_class, stage) %>%
    summarise(stage_avg = mean(count)) %>%
    ggplot(aes(x=stage, y=log1p(stage_avg), group = repeat_class)) +
    geom_line() +
    facet_wrap(.~ repeat_class, scales = "free") +
    theme_classic2(base_size = 20) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}




## FigS3A by SCAN-seq

if(F) {
  gene_list <- c("Tet1", "Tet2", "Tet3", "Tdg", "Dnmt1", "Uhrf1", "Dnmt3a", "Dnmt3b", "Dnmt3l")

  mouse_preim <- readRDS("data/processed/Mouse_preimplan.seurat.rds")

  stage_level <- c("M2","ZY","E2C", "L2C", "4C", "8C", "Morula", "ICM","TE")

  mouse_preim$Stage <- factor(mouse_preim$Stage, levels = stage_level)

  VlnPlot(mouse_preim, features = gene_list, group.by = "Stage")


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

  gene_list <- c("Tet1", "Tet2", "Tet3", "Tdg", "Dnmt1", "Uhrf1", "Dnmt3a", "Dnmt3b", "Dnmt3l")

  rna_mtx <- readRDS('data/processed/GSE136714.aggregated_gene_matrix.Rds')

  rna_mtx <- t(t(rna_mtx)/colSums(rna_mtx)*1e6)

  ### normalization (count -> FPKM) ###
  gene_anno <- readRDS('data/public/gencode.vM11.gene_length.Rds')

  ## check gene symbol & ensembl
  gene_sub <- gene_anno %>%
    filter(symbol %in% gene_list & symbol %in% rownames(rna_mtx)) %>%
    filter(!duplicated(symbol))

  gene_length <- gene_sub %>%
    dplyr::select(exonic.gene.sizes, symbol)
  rownames(gene_length) <- gene_length$symbol

  plot_mtx <- rna_mtx[gene_sub$symbol,]/gene_length[gene_sub$symbol,"exonic.gene.sizes"]*1e3

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

    # if(grepl("Tet", g)){
    if(T){
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
    summarise(mean=mean(value)) %>%
    tidyr::spread(stage, mean)
}

## FigS3A by Gao Shaorong (GSE98150)
if(F){

  gene_list <- c("Tet1", "Tet2", "Tet3", "Tdg", "Dnmt1", "Uhrf1", "Dnmt3a", "Dnmt3b", "Dnmt3l")

  gse_path <- "/mnt/f/Project/sc5hmC/data/ref/total_RNA/expr_mtx/mm10"
  rna_mtx <- read.table(file.path(gse_path, "GSE98150.totalRNA_FPKM.mm10.txt"), header = T)

  plot_df <- rna_mtx %>%
    filter(GeneSymbol %in% gene_list) %>%
    dplyr::select(-GeneLen) %>%
    reshape2::melt(id.vars = "GeneSymbol") %>%
    tidyr::separate(
      col = "variable",
      into = c("srr", "stage", NA),
      sep = "_") %>%
    filter(!grepl("6.5D|6.6D", stage))

  stage_level <- c("M2", "2C", "4C", "8C", "Morula", "ICM", "TE")

  plot_df$stage <- factor(plot_df$stage, levels = stage_level)

  todo_index <- 1:length(gene_list)
  barplot_fun <- function(x){
    g <- gene_list[x]
    p <- plot_df %>% filter(GeneSymbol==g) %>%
      ggbarplot(
        x = "stage", y = "value",fill="#6967aa",
        add = "mean_se", error.plot = "upper_errorbar")+
      ylim(min(plot_df$value), max(plot_df$value))+
      labs(y="Expression level", title=g)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.title = element_text(hjust = 0.5, face = "italic"))
    return(p)
  }

  p_list <- lapply(todo_index, barplot_fun)
  do.call("plot_grid", p_list)

  ggsave(
    "viz/FigS3A.5mC_related_gene_expression.GSR.pdf",
    width = 7.5, height = 5.7
  )

  plot_df %>%
    group_by(GeneSymbol, stage) %>%
    filter(grepl("Tet", GeneSymbol)) %>%
    summarise(mean=mean(value)) %>%
    tidyr::spread(stage, mean)
}

## FigS4D&E by Gao Shaorong (GSE98150)
if(F){

  gene_cor_res <- read.table(
    "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
    header = T, stringsAsFactors = F
  )

  pos_genes <- gene_cor_res %>%
    filter(pval < 0.05 & rho > 0) %>% pull(gene)
  neg_genes <- gene_cor_res %>%
    filter(pval < 0.05 & rho < 0) %>% pull(gene)

  gene_list <- c(pos_genes, neg_genes)

  gse_path <- "/mnt/f/Project/sc5hmC/data/ref/total_RNA/expr_mtx/mm10"
  rna_mtx <- read.table(file.path(gse_path, "GSE98150.totalRNA_FPKM.mm10.txt"), header = T)

  plot_df <- rna_mtx %>%
    filter(GeneSymbol %in% gene_list) %>%
    dplyr::select(-GeneLen) %>%
    reshape2::melt(id.vars = "GeneSymbol") %>%
    tidyr::separate(
      col = "variable",
      into = c("srr", "stage", NA),
      sep = "_") %>%
    filter(!grepl("6.5D|6.6D", stage)) %>%
    mutate(
      gene_group = case_when(
        GeneSymbol %in% pos_genes ~ "Positive",
        GeneSymbol %in% neg_genes ~ "Negative",
        TRUE ~ "Others"
      )
    )

  stage_level <- c("M2", "2C", "4C", "8C", "Morula", "ICM", "TE")

  plot_df$stage <- factor(plot_df$stage, levels = stage_level)

  too_high <- quantile(plot_df$value, probs = 0.90)

  pal <- c("Positive"="#e51d17", "Negative"="#3d51a2")

  plot_df %>%
    filter(value < too_high) %>%
    ggplot(aes(x=stage, y=value, color=gene_group)) +
    geom_violin(scale = "width") +
    facet_wrap(gene_group~., scales = "free_y") +
    scale_color_manual(values = pal) +
    theme_classic2(base_size = 20) +
    theme(aspect.ratio = 1,
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  ggsave(
    "viz/FigS4D&E.5mC_related_gene_expression.GSR.240510.pdf",
    width = 7.5, height = 5.7
  )

}

## FigS4D&E by SCAN-seq
if(F){

  gene_cor_res <- read.table(
    "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
    header = T, stringsAsFactors = F
  )

  pos_genes <- gene_cor_res %>%
    filter(pval < 0.05 & rho > 0) %>% pull(gene)
  neg_genes <- gene_cor_res %>%
    filter(pval < 0.05 & rho < 0) %>% pull(gene)

  gene_list <- c(pos_genes, neg_genes)

  mouse_preim <- readRDS("data/processed/Mouse_preimplan.seurat.rds")

  mouse_preim$Stage_new <- case_when(
    mouse_preim$Stage %in% c("E2C", "L2C") ~ "2C",
    mouse_preim$Stage %in% c("ICM", "TE") ~ "Blast",
    TRUE ~ mouse_preim$Stage) %>%
    setNames(names(mouse_preim$Stage))

  stage_level <- c("M2","ZY","2C", "4C", "8C", "Morula", "Blast")

  mouse_preim$Stage_new <- factor(
    mouse_preim$Stage_new, levels = stage_level
  )
  table(mouse_preim$Stage_new)

  # mouse_preim <- AddModuleScore(
  #   mouse_preim,
  #   features = list("pos"=pos_genes, "neg"=neg_genes),
  #   name = c("pos_cor", "neg_cor")
  # )
  #
  # VlnPlot(mouse_preim, features = c("pos_cor1", "neg_cor2"), group.by = "Stage")

  plot_df <- mouse_preim@assays$RNA@data %>%
    as.matrix() %>% .[intersect(gene_list, rownames(mouse_preim)), ] %>%
    reshape2::melt() %>%
    left_join(
      mouse_preim@meta.data[,c("Cell_id", "Stage_new")],
      by=c("Var2"="Cell_id")) %>%
    mutate(
      gene_group = case_when(
        Var1 %in% pos_genes ~ "Positive",
        Var1 %in% neg_genes ~ "Negative",
        TRUE ~ "Others"
      )
    )
  plot_df$Stage_new <- factor(plot_df$Stage_new, levels = stage_level)

  plot_df %>%
    filter(value < 1) %>%
    ggplot(aes(x=Stage_new, y=value, color=gene_group)) +
    geom_violin(scale = "width") +
    facet_wrap(gene_group~.) +
    scale_color_manual(values = pal) +
    theme_classic2(base_size = 20) +
    theme(aspect.ratio = 1,
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}



