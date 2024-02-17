## library packages
library(data.table)
library(dplyr)
library(ggplot2)
library(Biobase)
## input
# setwd("/mnt/f/Project/sc5hmC/data/ref/5mC_5hmC_genebody/")
# mc_mtx <- read.table("mC_gene_ratio.tsv", header = T, stringsAsFactors = F)
# hmc_mtx <- read.table("1111_hmC_gene_ratio.tsv", header = T, stringsAsFactors = F)
# common_col <- c("chr"="chr","start"="start","end"="end","gene"="gene")
setwd("/mnt/f/Project/sc5hmC/data/")

hmc_mtx <- lapply(
  list.files("DHMR/hmC_genebody/", full.names = T),
  function(x){
    df <- read.table(x, header = F, stringsAsFactors = F)
    sp <- strsplit(basename(x), "[.]")[[1]][1]
    # sp <- strsplit(sp, "_")[[1]][2]
    colnames(df) <- c("chr", "start", "end", "gene", sp)
    df <- df %>% .[,c("gene", sp)]
    print(sapply(df, class))
    return(df)
  }
)
hmc_mtx <- hmc_mtx %>% purrr::reduce(function(x, y)
  base::merge(x, y, by = c("gene"))) %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.))
head(hmc_mtx)

# range(hmc_mtx)

##--------------------5mC--------------------
mc_mtx_lj <-
  read.table(
    "ref/5mC_5hmC_genebody/mC_gene_ratio.tsv",
    header = T,
    stringsAsFactors = F
  )
mc_mtx_lj <- mc_mtx_lj %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.)) %>%
  dplyr::select(gene, oocyte, sperm)
# mc_mtx <- mc_mtx %>%
#   filter(!duplicated(gene)) %>%
#   tibble::rownames_to_column("RowNames") %>%
#   select(-chr, -start, -end, -RowNames) %>%
#   tibble::column_to_rownames("gene") %>%
#   t %>% scale() %>% t

mc_mtx_gsr <- lapply(
  list.files("DHMR/mC_genebody/", full.names = T),
  function(x){
    df <- read.table(x, header = F, stringsAsFactors = F)
    sp <- strsplit(basename(x), "[.]")[[1]][1]
    sp <- strsplit(sp, "_")[[1]][2]
    colnames(df) <- c("chr", "start", "end", "gene", sp)
    df <- df %>% .[,c("gene", sp)]
    print(sapply(df, class))
    return(df)
  }
)

mc_mtx_gsr <-
  mc_mtx_gsr %>% purrr::reduce(function(x, y)
    base::merge(x, y, by = c("gene"))) %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.))

mc_mtx <- inner_join(
  mc_mtx_lj, mc_mtx_gsr,by = "gene")
head(mc_mtx)

##--------------method 1: geom_smooth---------------
plot_df <- list(mc_mtx, hmc_mtx) %>%
  purrr::reduce(function(x,y)inner_join(x,y,by=common_col, suffix = c("-mC", "-hmC")))

# coeff1 <- max(plot_df$X4c)/max(plot_df$`4C`)
# coeff2 <- max(plot_df$four_cell)/max(plot_df$`4C`)

stage_list <- colnames(plot_df)[grepl("-", colnames(plot_df))]
stage_list <- sapply(strsplit(stage_list, "-"), "[[", 1) %>% unique()

todo_index <- 1:length(stage_list)
plot_5mC_5hmC <- function(x){
  stg <- stage_list[x]
  coeff3 <- max(plot_df[,paste0(stg, "-mC")], na.rm = T)/
    max(plot_df[,paste0(stg, "-hmC")], na.rm = T)
  p <- plot_df %>%
    arrange(desc(get(paste0(stg, "-hmC")))) %>%
    mutate(rank=1:nrow(plot_df)) %>%
    ggplot(aes(x=rank))+
    ## red
    geom_smooth(aes(y=get(paste0(stg, "-hmC"))), color="#ff0000", size=2.5)+
    ## purple
    geom_smooth(aes(y=get(paste0(stg, "-mC"))/coeff3), color="#a020f0", se = F, size=2.5) +
    ## green
    # geom_smooth(aes(x=rank,y=four_cell/coeff2), color="green", size=2.5)+
    scale_y_continuous(
      name = "5hmC/C",
      sec.axis = sec_axis(trans = ~. *coeff3, name = "5mC/C"))+
    labs(x="Ranked 5hmC levels", title = stg)+
    theme_bw(base_size = 20)+
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = .5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  return(p)
}

p_list <- lapply(todo_index, plot_5mC_5hmC)

library(cowplot)
p <- plot_df %>% select(`two_cell-hmC`, `two_cell-hmC`) %>%
  reshape2::melt() %>% head() %>% mutate(group=gsub("two_cell-", "5", variable)) %>%
  ggplot(aes(x=group, y=value, color=group))+geom_line(size=2.5)+
  scale_color_manual(values = c("5hmC"="#ff0000", "5mC"="#a020f0"))
legend <- get_legend(p)
p_list[[length(p_list)+1]] <- legend

do.call("plot_grid", p_list)
ggsave("F:/Project/sc5hmC/manuscript/Fig3A.pdf", width = 15, height = 8)

##--------------method 2: bin the genes by expression level--------------
library(Seurat)
library(egg)
mouse_preim <- readRDS("ref/scan-seq/data/Gene/Mouse_preimplan.seurat.rds")
rna_mtx <- (mouse_preim@assays$RNA@data %>% as.matrix() %>% exp)-1

rna_mtx <- rna_mtx %>%
  reshape2::melt() %>%
  left_join(
    mouse_preim@meta.data %>%
      select(c("Cell_id", "Stage")),
    by = c("Var2" = "Cell_id")) %>%
  mutate(
    Stage = case_when(
      Stage %in% c("E2C", "L2C") ~ "2C",
      Stage %in% c("ICM", "TE") ~ "Blast",
      TRUE ~ Stage )) %>%
  dplyr::group_by(Stage, Var1) %>%
  # dplyr::summarise(mean_expr=mean(value, na.rm=T)) %>%
  dplyr::summarise(mean_expr=log1p(mean(value, na.rm=T))) %>%
  tidyr::spread(key = Stage, value = mean_expr)
#
# rna_mtx_com$`2C` <- rowMeans(rna_mtx_com[,c("E2C", "L2C")])
# rna_mtx_com$`Blast` <- rowMeans(rna_mtx_com[,c("TE", "ICM")])
rna_mtx <- rna_mtx %>% dplyr::select(Var1,`M2`, `ZY`,`2C`,`4C`, `8C`, `Morula`, `Blast`)
colnames(rna_mtx) <- c("gene","Oocyte_rna", "Zygote_rna","2C_rna", "4C_rna","8C_rna", "Morula_rna", "Blast_rna")
# rna_mtx <- rna_mtx %>% tibble::column_to_rownames("gene")
head(rna_mtx)

# mouse_preim <- readRDS("/mnt/f/Project/sc5hmC/data/ref/scan-seq/data/Gene/Mouse_preimplan.seurat.rds")
# rna_mtx <- (mouse_preim@assays$RNA@data %>% as.matrix() %>% exp)-1
#
# common_genes <- list(rownames(rna_mtx), mc_mtx$gene, hmc_mtx$gene) %>%
#   purrr::reduce(function(x,y)intersect(x,y))
#
# rna_mtx_com <- rna_mtx[common_genes,] %>% reshape2::melt() %>%
#   left_join(mouse_preim@meta.data[,c("Cell_id", "Stage")], by=c("Var2"="Cell_id"))%>%
#   dplyr::group_by(Stage, Var1) %>%
#   # dplyr::summarise(mean_expr=mean(value, na.rm=T)) %>%
#   dplyr::summarise(mean_expr=log1p(mean(value, na.rm=T))) %>%
#   tidyr::spread(key = Stage, value = mean_expr)
# rna_mtx_com$`2C` <- rowMeans(rna_mtx_com[,c("E2C", "L2C")])
# rna_mtx_com$`Blast` <- rowMeans(rna_mtx_com[,c("TE", "ICM")])
# rna_mtx_com <- rna_mtx_com %>% dplyr::select(Var1, `2C`,`4C`,`M2`, `Blast`)
# colnames(rna_mtx_com) <- c("gene", "two_cell-rna", "four_cell-rna","oocyte-rna", "blast-rna")

# mc_mtx_com <- mc_mtx %>% filter(gene %in% common_genes)
# hmc_mtx_com <- hmc_mtx %>% filter(gene %in% common_genes)
head(mc_mtx)
head(hmc_mtx)

colnames(mc_mtx) <- c("gene", "Oocyte", "Sperm", "Zygote", "2C", "4C", "8C", "Morula", "Blast")

plot_df <- inner_join(
  mc_mtx, hmc_mtx, by=c("gene"), suffix = c("_mC", "_hmC")) %>%
  inner_join(rna_mtx, by=c("gene"="gene"))

# plot_df <- list(rna_mtx_com, mc_mtx_com, hmc_mtx_com) %>%
#   purrr::reduce(function(x,y)inner_join(x,y,by=c("gene"="gene")))

stage_list <- colnames(plot_df)[grepl("_", colnames(plot_df))]
stage_list <- sapply(strsplit(stage_list, "_"), "[[", 1) %>% sort()
stage_list <- names(table(stage_list))[table(stage_list)==3]
stage_list

library(scales)
library(forecast)
todo_index <- 1:length(stage_list)
plot_5mC_5hmC_RNA <- function(x){
  # stg <- stage_list[x]
  ## triple-omics
  # plot_df %>%
  #   mutate(gene_bin = cut(get(paste0(stg, "-rna")), breaks=300)) %>%
  #   group_by(gene_bin) %>%
  #   summarise(
  #     mean_rna=mean(get(paste0(stg, "-rna")), na.rm = T),
  #     mean_hmC=mean(get(paste0(stg, "-hmC")), na.rm = T),
  #     mean_mC=mean(get(paste0(stg, "-mC")), na.rm = T)) %>%
  #   arrange(desc(gene_bin)) %>%
  #   mutate(rank=order(gene_bin, decreasing = T)) %>%
  #   mutate(coeff1=max(mean_mC)/max(mean_hmC),
  #          coeff2=max(mean_rna)/max(mean_hmC)) %>%
  #   ggplot(aes(x=rank))+
  #   geom_smooth(aes(y=mean_hmC), color="#a020f0", size=2.5, se=F) + ## purple
  #   geom_smooth(aes(y=mean_mC/coeff1), color="green", size=2.5, se=F)+ ## green
  #     geom_point(aes(y=mean_hmC), color="#a020f0", size=2.5) + ## purple
  #     geom_point(aes(y=mean_mC/coeff1), color="green", size=2.5)+ ## green
  #   # geom_point(aes(y=mean_rna/coeff2), color="#ff0000", size=2.5)+ ## red
  #   # scale_y_continuous(
  #   #   name = "5hmC or 5mC level",
  #   #   sec.axis = sec_axis(trans = ~. /coeff2, name = "RNA expression level"))+
  #   labs(x="Ranked genes by expression level", title = stg)+
  #   theme_bw(base_size = 20)+
  #   theme(panel.grid = element_blank(),
  #         plot.title = element_text(hjust = .5),
  #         axis.ticks.x = element_blank(),
  #         axis.text.x = element_blank())
  #
  # ## RNA expression as heatmap
  # tri_omics_df <- plot_df %>%
  #   mutate(gene_bin = cut(get(paste0(stg, "-rna")), breaks=300)) %>%
  #   group_by(gene_bin) %>%
  #   summarise(
  #     mean_rna=mean(get(paste0(stg, "-rna")), na.rm = T),
  #     mean_hmC=mean(get(paste0(stg, "-hmC")), na.rm = T),
  #     mean_mC=mean(get(paste0(stg, "-mC")), na.rm = T)) %>%
  #   arrange(desc(gene_bin)) %>%
  #   mutate(rank=order(gene_bin, decreasing = T)) %>%
  #   mutate(coeff1=max(mean_mC)/max(mean_hmC))
  #
  # cor_res_1 <- cor.test(tri_omics_df$mean_hmC, tri_omics_df$mean_mC, method = "spearman")
  # rho1 <- round(cor_res$estimate,2)
  # p1 <- format(cor_res_1$p.value, digits=3)
  #
  # cor_res_2 <- cor.test(tri_omics_df$mean_hmC, tri_omics_df$mean_rna, method = "spearman")
  # rho2 <- round(cor_res_2$estimate,2)
  # p2 <- format(cor_res_2$p.value, digits=3)
  #
  # cor_res_3 <- cor.test(tri_omics_df$mean_mC, tri_omics_df$mean_rna, method = "spearman")
  # rho3 <- round(cor_res_3$estimate,2)
  # p3 <- format(cor_res_3$p.value, digits=3)
  #
  # ratio <- diff(range(tri_omics_df$rank))/diff(range(tri_omics_df$mean_mC))
  #
  # title <- paste0(
  #   "Stage:", stg, "\n",
  #   "5hmC vs 5mC:rho=",  rho1, ", pvalue=", p1, "\n",
  #   "5hmC vs RNA:rho=",  rho2, ", pvalue=", p2,"\n",
  #   "5mC vs RNA:rho=",  rho3, ", pvalue=", p3
  # )
  #
  # p1 <- tri_omics_df %>%
  #   ggplot(aes(x=rank))+
  #   geom_smooth(aes(y=mean_hmC), color="#a020f0", size=2.5, se=F,span=2) + ## purple
  #   geom_smooth(aes(y=mean_mC/coeff1), color="green", size=2.5, se=F,span=2)+ ## green
  #   # geom_point(aes(y=mean_rna/coeff2), color="#ff0000", size=2.5)+ ## red
  #   scale_y_continuous(
  #     name = "5hmC or 5mC level",
  #     sec.axis = sec_axis(
  #       trans = ~. *unique(tri_omics_df$coeff1),
  #       name = "RNA expression level"))+
  #   scale_x_continuous(expand = c(0,0))+
  #   # coord_fixed()+
  #   # labs(title = paste0(stg, " rho:", rho))+
  #   labs(title = title)+
  #   theme_bw(base_size = 15)+
  #   theme(panel.grid = element_blank(),
  #         plot.title = element_text(hjust = 0, size=10),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.text.x = element_blank())
  #
  # p2 <- tri_omics_df %>%
  #   ggplot(aes(x=rank, y="rna", fill=mean_rna))+
  #   geom_tile()+labs(y="")+
  #   scale_x_continuous(expand = c(0,0))+
  #   labs(x="Ranked genes by expression level")+
  #   # scale_fill_viridis_c()+
  #   scale_fill_gradient2(low = "#075AFF",mid = "#FFFFCC",high = "#FF0000")+
  #   theme_minimal(base_size = 15)+
  #   theme(panel.grid = element_blank(),
  #         legend.position = "bottom",
  #         axis.title.y = element_blank(),
  #         axis.text = element_blank())
  #
  # # plot_grid(p1,p2, nrow = 2,rel_heights = c(5,2))
  # pl <- ggarrange(p1,p2,ncol = 1,heights = c(5,1))
  stg <- stage_list[x]

  index_mc <- order(plot_df[,paste0(stg, "_mC")])
  index_hmc <- order(plot_df[,paste0(stg, "_hmC")])
  index_rna <- order(plot_df[,paste0(stg, "_rna")])

  ## moving average
  tri_omics_df <- plot_df %>%
    # mutate(gene_bin = cut(get(paste0(stg, "-rna")), breaks=300)) %>%
    # group_by(gene_bin) %>%
    # summarise(
    #   mean_rna=mean(get(paste0(stg, "-rna")), na.rm = T),
    #   mean_hmC=mean(get(paste0(stg, "-hmC")), na.rm = T),
    #   mean_mC=mean(get(paste0(stg, "-mC")), na.rm = T)) %>%
    # dplyr::arrange(get(paste0(stg, "-rna"))) %>%
    dplyr::mutate(
      mean_mC=ma(plot_df[index_mc,paste0(stg, "_mC")], order = 200, centre = F),
      # mean_mC=ma(plot_df[index_mc,paste0(stg, "_mC")], order = 10, centre = F),
      mean_hmC=ma(plot_df[index_hmc,paste0(stg, "_hmC")], order = 200, centre = F),
      # mean_hmC=ma(plot_df[index_hmc,paste0(stg, "_hmC")], order = 10, centre = F),
      mean_rna=ma(plot_df[index_rna,paste0(stg, "_rna")], order = 10, centre = F)
    ) %>%
    na.omit() %>%
    dplyr::arrange(desc(mean_rna)) %>%
    dplyr::mutate(rank=order(mean_rna, decreasing = T)) %>%
    dplyr::mutate(coeff1=max(mean_mC)/max(mean_hmC))

  cor_res_1 <- cor.test(plot_df[,paste0(stg, "_hmC")],plot_df[,paste0(stg, "_mC")], method = "spearman")
  rho1 <- round(cor_res_1$estimate,2)
  p1 <- format(cor_res_1$p.value, digits=3)

  cor_res_2 <- cor.test(plot_df[,paste0(stg, "_hmC")],plot_df[,paste0(stg, "_rna")], method = "spearman")
  rho2 <- round(cor_res_2$estimate,2)
  p2 <- format(cor_res_2$p.value, digits=3)

  cor_res_3 <- cor.test(plot_df[,paste0(stg, "_mC")],plot_df[,paste0(stg, "_rna")], method = "spearman")
  rho3 <- round(cor_res_3$estimate,2)
  p3 <- format(cor_res_3$p.value, digits=3)

  ratio <- diff(range(tri_omics_df$rank))/diff(range(tri_omics_df$mean_mC))

  title <- paste0(
    "Stage:", stg, "\n",
    "5hmC vs 5mC:rho=",  rho1, ", pvalue=", p1, "\n",
    "5hmC vs RNA:rho=",  rho2, ", pvalue=", p2,"\n",
    "5mC vs RNA:rho=",  rho3, ", pvalue=", p3
  )

  p1 <- tri_omics_df %>%
    ggplot(aes(x=rank))+
    geom_smooth(aes(y=mean_hmC), color="#a020f0", size=2.5, se=F,span=1) + ## purple
    geom_smooth(aes(y=mean_mC/coeff1), color="#58b050", size=2.5, se=F,span=1)+ ## green
    # geom_point(aes(y=mean_hmC), color="#a020f0", size=1, alpha=.5) + ## purple
    # geom_point(aes(y=mean_mC/coeff1), color="green", size=1, alpha=.5)+ ## green
    # geom_point(aes(y=mean_rna/coeff2), color="#ff0000", size=2.5)+ ## red
    scale_y_continuous(
      name = "5hmCpG level",
      sec.axis = sec_axis(
        trans = ~. *unique(tri_omics_df$coeff1),
        name = "5mCpG level"))+
    scale_x_continuous(expand = c(0,0))+
    # coord_fixed()+
    # labs(title = paste0(stg, " rho:", rho))+
    labs(title = title)+
    theme_bw(base_size = 15)+
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0, size=10),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())

  p2 <- tri_omics_df %>%
    ggplot(aes(x=rank, y="rna", fill=mean_rna))+
    geom_tile()+labs(y="")+
    scale_x_continuous(expand = c(0,0))+
    labs(x="Ranked genes by expression level")+
    # scale_fill_viridis_c()+
    scale_fill_gradient2(low = "#075AFF",mid = "#FFFFCC",high = "#FF0000")+
    theme_minimal(base_size = 15)+
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          axis.title.y = element_blank(),
          axis.text = element_blank())

  # plot_grid(p1,p2, nrow = 2,rel_heights = c(5,2))
  # pl <- ggarrange(p1,p2,ncol = 1,heights = c(5,1))
  pl <- cowplot::plot_grid(p1,p2,ncol = 1, rel_heights = c(0.8,0.2))
  # pl
  return(pl)
}

p_list <- lapply(todo_index, plot_5mC_5hmC_RNA)

p <- plot_df %>% select(`2C_hmC`, `2C_hmC`) %>%
  reshape2::melt() %>% mutate(group=gsub("2C_", "5", variable)) %>%
  ggplot(aes(x=group, y=value, color=group))+geom_line(size=2.5)+
  # scale_color_manual(values = c("5hmC"="#a020f0", "5mC"="green"))
  scale_color_manual(values = c("5hmC"="#a020f0", "5mC"="#58b050"))
legend <- cowplot::get_legend(p)
p_list[[length(p_list)+1]] <- legend

cowplot::plot_grid(plotlist = p_list, nrow = 2)
# do.call("plot_grid", c(p_list,nrow=1))

ggsave("Fig4A.add_RNA.span_2.1224.pdf", width = 40, height = 5)

## correlation in the pos or neg genes

gene_cor_res <- read.table("Gene_wise_cor.ZY_Bla.5hmC_5mC.tsv", header = T, stringsAsFactors = F)

pos_genes <- gene_cor_res %>% filter(pval < 0.05 & rho > 0) %>% pull(gene)
neg_genes <- gene_cor_res %>% filter(pval < 0.05 & rho < 0) %>% pull(gene)

plot_orig_df <- plot_df

plot_df <- plot_orig_df %>% filter(gene %in% neg_genes)

p_list <- lapply(todo_index, plot_5mC_5hmC_RNA)

p <- plot_df %>% select(`2C_hmC`, `2C_hmC`) %>%
  reshape2::melt() %>% mutate(group=gsub("2C_", "5", variable)) %>%
  ggplot(aes(x=group, y=value, color=group))+geom_line(size=2.5)+
  # scale_color_manual(values = c("5hmC"="#a020f0", "5mC"="green"))
  scale_color_manual(values = c("5hmC"="#a020f0", "5mC"="#58b050"))
legend <- cowplot::get_legend(p)
p_list[[length(p_list)+1]] <- legend

cowplot::plot_grid(plotlist = p_list, nrow = 2)

p1 <- plot_df %>%
  ggplot(aes(x=`2C_rna`, y=`2C_hmC`))+
  geom_point(alpha=.3, color="grey") +
  geom_density_2d(color="black") +
  ylim(c(0,0.06)) +xlim(c(0,4))  +
  labs(title = "all the genes") +
  theme_bw() +
  theme(aspect.ratio = 1)

p1 <- ggExtra::ggMarginal(p1, type = "histogram", fill="grey")

p2 <- plot_df %>%
  filter(gene %in% pos_genes) %>%
  ggplot(aes(x=`2C_rna`, y=`2C_hmC`))+
  geom_point(alpha=.3, color="#df6c77") +
  geom_density_2d(color="black") +
  xlim(c(0,3)) + ylim(c(0,0.05)) +
  labs(title = "postively correlated genes") +
  theme_bw() +
  theme(aspect.ratio = 1)

p2 <- ggExtra::ggMarginal(p2, type = "histogram", fill="#df6c77")

p3 <- plot_df %>%
  filter(gene %in% neg_genes) %>%
  ggplot(aes(x=`2C_rna`, y=`2C_hmC`))+
  geom_point(alpha=.3, color="#80B1D3") +
  geom_density_2d(color="black") +
  xlim(c(0,3)) + ylim(c(0,0.05)) +
  labs(title = "negatively correlated genes") +
  theme_bw() +
  theme(aspect.ratio = 1)

p3 <- ggExtra::ggMarginal(p3, type = "histogram", fill="#80B1D3")

cowplot::plot_grid(p1,p2,p3, ncol = 1)

########################################
##   correlation in the genome-wide   ##
########################################
##--------------packages---------------
library(data.table)
library(dplyr)
library(ggplot2)
##--------------input---------------
## merge data of early 2 cell and late 2 cell

setwd("F:/Project/sc5hmC/data/ref/WBGS_GSR/")
merge_fun <- function(dir1, dir2, out_name){
  input_file1 <- fread(dir1)
  input_file2 <- fread(dir2)

  out_file <- inner_join(
    out_file1, out_file2,by=c("V1"="V1", "V2"="V2", "V3"="V3")
  )
  out_file$cov <- as.integer(out_file$V4.x+out_file$V4.y)
  out_file$meth <- as.integer(out_file$V5.x+out_file$V5.y)
  out_file <- out_file[,c(1:3, 8,9)]

  out_dir <- strsplit(out_name,"/")[[1]][1]
  dir.create(out_dir)
  write.table(out_file, file = out_name,col.names = F, row.names = F, quote = F, sep = "\t")
}

merge_fun(
  dir1 = "other/GSM2588709_2cellEarly/GSM2588709_2cellEarly.1X50.aim.1000_1X.bed",
  dir2 = "other/GSM2588710_2cellLate/GSM2588710_2cellLate.1X50.aim.1000_1X.bed",
  out_name = "GSM2588709.5_2cell/GSM2588709.5_2cell.1X50.aim.1000_1X.bed"
)

merge_fun(
  dir1 = "other/GSM2588714_ICM/GSM2588714_ICM.1X50.aim.1000_1X.bed",
  dir2 = "other/GSM2588715_TE/GSM2588715_TE.1X50.aim.1000_1X.bed",
  out_name = "GSM2588714.5_Blast/GSM2588714.5_Blast.1X50.aim.1000_1X.bed"
)

##--------------compare 5hmC data and 5mC data tile-level-------------

stage_list <- c("2cell", "4cell", "")

mc_dt <- fread("GSM2588709.5_2cell/GSM2588709.5_2cell.1X50.aim.1000_1X.bed", stringsAsFactors = F)
hmc_dt <- fread("../merged_5hmC/2C/2C.1X50.aim.1000_1X.bed", stringsAsFactors = F)

mc_dt <- mc_dt %>% mutate(mc_ratio=V5/V4) %>% select(V1,V2,V3,mc_ratio)
hmc_dt <- hmc_dt %>% mutate(hmc_ratio=V5/V4) %>% select(V1,V2,V3,hmc_ratio)
merge_dt <- inner_join(mc_dt, hmc_dt, by=c("V1"="V1", "V2"="V2", "V3"="V3"))
merge_dt<- merge_dt %>% arrange(desc(hmc_ratio))

merge_dt[,] %>%
  filter(hmc_ratio <= 0.02) %>%
  filter(!duplicated(hmc_ratio)) %>%
  arrange(desc(hmc_ratio)) %>%
  # filter(mc_ratio <= 0.5 & mc_ratio > 0.3) %>%
  mutate(rank=order(hmc_ratio, decreasing = T)) %>%
  ggplot()+
  ## red
  geom_smooth(aes(x=rank,y=hmc_ratio), color="#ff0000", size=2.5)+
  ## purple
  geom_smooth(aes(x=rank, y=mc_ratio/50), color="#a020f0", se = F, size=2.5) +
  scale_y_continuous(
    name = "5hmC/C",
    sec.axis = sec_axis(trans = ~. *50, name = "5mC/C"))+
  labs(x="Ranked 5hmC levels", title = "2C")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = .5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())



if (F) {
  input_tmp <- plot_df %>% filter(`blast-rna` <=as.numeric(quantile(plot_df$`blast-rna`, probs = c(0.9))))
  glmm<-glm(formula = `blast-mC`~`blast-rna`,data = input_tmp, family = Gamma)
  plot(glmm)

  input_tmp <- plot_df %>% filter(`blast-rna` <=as.numeric(quantile(plot_df$`blast-rna`, probs = c(0.9))))
  m <- glm(formula = `blast-mC`~`blast-rna`,data = input_tmp, family = nbiom)
  plot(m)

  # input_tmp <- plot_df %>% filter(`blast-hmC` > 0)
  input_tmp <- plot_df %>% mutate(gene_len=end-start)

  plot(input_tmp$`blast-hmC`, input_tmp$`blast-mC`, cex=.1)

  m <- lm(`blast-hmC` ~ `blast-mC` + log1p(gene_len) , data = input_tmp)
  # m <- MASS::glm.nb(`blast-hmC` ~ `blast-mC`, data = input_tmp)
  hist(residuals(m),breaks = 100)

  anova(m)
  summary(m)

  input_tmp <- plot_df %>% filter(`blast-rna` <=as.numeric(quantile(plot_df$`blast-rna`, probs = c(0.9))))
  lmm<-lm(formula = `blast-mC`~`blast-rna`,data = input_tmp)
  plot(lmm)

  lmm<-lm(formula = `blast-rna`~`blast-mC`+`blast-hmC`,data = input_tmp)
  anova(lmm)
  summary(lmm)

  lmm<-lm(formula = `blast-rna`~`blast-hmC`+`blast-mC`,data = input_tmp)
  summary(lmm)

  library(forecast)

  index<-order(plot_df$`blast-rna`)
  plot(plot_df$`blast-rna`[index] %>% log1p, plot_df$`blast-hmC`[index])
  lines(x=plot_df$`blast-rna`[index] %>% log1p,
        y=ma(plot_df$`blast-hmC`[index]%>% log1p,order=500) ,
        cex=6,
        col='red')

  plot(plot_df$`blast-rna`[index] %>% log1p, plot_df$`blast-mC`[index])
  lines(x=plot_df$`blast-rna`[index] %>% log1p,
        y=ma(plot_df$`blast-mC`[index]%>% log1p,order=500) ,
        cex=6,
        col='red')

  cor.test(plot_df$`blast-rna`,plot_df$`blast-hmC`,method = 'spearman')

  cor.test(plot_df$`blast-rna`[index] %>% log1p,
           ma(plot_df$`blast-hmC`[index]%>% log1p,order=500))

  cor.test(plot_df$`blast-rna`[index] %>% log1p,
           ma(plot_df$`blast-mC`[index]%>% log1p,order=500))

  index<-order(input_tmp$`blast-rna`)
  plot(input_tmp$`blast-rna`[index] %>% log1p, input_tmp$`blast-hmC`[index])
  lines(x=input_tmp$`blast-rna`[index] %>% log1p,
        y=ma(input_tmp$`blast-hmC`[index]%>% log1p,order=500) ,
        cex=6,
        col='red')

  blast.mc.res <- lm(formula =`blast-mC` ~ `blast-hmC`,data = plot_df)$residual
  blast.hmc.res <- lm(formula =`blast-hmC` ~ `blast-mC` ,data = plot_df)$residual
  plot(blast.hmc.res,plot_df$`blast-rna` %>% log1p)

  # 5mc residual ~ rna
  index<-order(plot_df$`blast-rna`)
  plot(plot_df$`blast-rna`[index] %>% log1p, blast.mc.res[index])
  lines(x=plot_df$`blast-rna`[index] %>% log1p,
        y=ma(blast.mc.res[index]%>% log1p,order=100) ,
        cex=6,
        col='red')
  cor.test(blast.mc.res,plot_df$`blast-rna`,method = 'spearman')
  #cor.test(plot_df$`blast-mC`-plot_df$`blast-hmC`,plot_df$`blast-rna`,method = 'spearman')
  cor.test(plot_df$`blast-hmC`,plot_df$`blast-rna`,method = 'spearman')

  plot(plot_df$`blast-rna`[index] %>% log1p, blast.hmc.res[index])
  lines(x=plot_df$`blast-rna`[index] %>% log1p,
        y=ma(blast.hmc.res[index]%>% log1p,order=100) ,
        cex=6,
        col='red')

  cor.test(plot_df$`blast-rna` %>% log1p, plot_df$`blast-hmC`,method = 'spearman')

  cor.test(plot_df$`blast-rna` %>% log1p, plot_df$`blast-mC` - plot_df$`blast-hmC`,method = 'spearman')
  lm(log1p(`blast-rna`) ~ `blast-mC` - `blast-hmC`, plot_df) %>% summary
  lm(log1p(`blast-rna`) ~ `blast-hmC` - `blast-mC`, plot_df) %>% summary
  lm(log1p(`blast-rna`) ~ `blast-mC`, plot_df) %>% summary
  cor.test(log1p(plot_df$`blast-rna`),plot_df$`blast-mC`,method = 'spearman')

  plot(plot_df$`blast-mC`,plot_df$`blast-hmC`)
  abline(a=1,b=0)
}
