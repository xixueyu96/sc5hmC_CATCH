library(scales)
library(forecast)
library(dplyr)
library(ggplot2)

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

stage_list <- colnames(plot_df)[grepl("_", colnames(plot_df))]
stage_list <- sapply(strsplit(stage_list, "_"), "[[", 1) %>% sort()
stage_list <- names(table(stage_list))[table(stage_list)==3]
stage_list

plot_multiomic_cor <- function(x){

  # stg <- stage_list[x]
  stg <- x

  index_mc <- order(plot_df[,paste0(stg, "_mC")])
  index_hmc <- order(plot_df[,paste0(stg, "_hmC")])
  index_rna <- order(plot_df[,paste0(stg, "_rna")])

  ## moving average
  tri_omics_df <- plot_df %>%
    mutate(
      mean_mC=ma(plot_df[index_mc,paste0(stg, "_mC")], order = 200, centre = F),
      # mean_mC=ma(plot_df[index_mc,paste0(stg, "_mC")], order = 10, centre = F),
      mean_hmC=ma(plot_df[index_hmc,paste0(stg, "_hmC")], order = 200, centre = F),
      # mean_hmC=ma(plot_df[index_hmc,paste0(stg, "_hmC")], order = 10, centre = F),
      mean_rna=ma(plot_df[index_rna,paste0(stg, "_rna")], order = 10, centre = F)
    ) %>% na.omit() %>%
    arrange(desc(mean_rna)) %>%
    mutate(rank=order(mean_rna, decreasing = T)) %>%
    mutate(coeff1=max(mean_mC)/max(mean_hmC))

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
    scale_y_continuous(
      name = "5hmCpG level",
      expand = c(0,0),
      sec.axis = sec_axis(
        trans = ~. *unique(tri_omics_df$coeff1),
        name = "5mCpG level"))+
    scale_x_continuous(expand = c(0,0))+
    labs(title = title)+
    theme_bw(base_size = 15)+
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          plot.title = element_text(hjust = 0, size=10),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())

  p2 <- tri_omics_df %>%
    ggplot(aes(x=rank, y="rna", fill=mean_rna))+
    geom_tile()+labs(y="")+
    scale_x_continuous(expand = c(0,0))+
    labs(x="Descending ordered genes")+
    # scale_fill_viridis_c()+
    scale_fill_gradient2(low = "#075AFF",mid = "#FFFFCC",high = "#FF0000")+
    theme_minimal(base_size = 15)+
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          axis.title.y = element_blank(),
          axis.text = element_blank())

  # plot_grid(p1,p2, nrow = 2,rel_heights = c(5,2))
  pl <- egg::ggarrange(p1, p2, ncol = 1, heights = c(0.8,0.2))
  # pl <- cowplot::plot_grid(p1,p2,ncol = 1, rel_heights = c(0.8,0.2), align = "hv")
  # pl
  return(pl)
}

# plot_multiomic_cor("Zygote") # test

p_list <- lapply(
  c("Oocyte",
    "Zygote",
    "Two.cell",
    "Four.cell",
    "Eight.cell",
    "Morula",
    "Blast"),
  plot_multiomic_cor
)

p <- plot_df %>%
  dplyr::select(Zygote_mC, Zygote_hmC) %>%
  reshape2::melt() %>%
  mutate(group=gsub("Zygote_", "5", variable)) %>%
  ggplot(aes(x=group, y=value, color=group)) +
  geom_line(size=2.5)+
  scale_color_manual(
    values = c("5hmC"="#a020f0", "5mC"="#58b050")
  )
legend <- cowplot::get_legend(p)
p_list[[length(p_list)+1]] <- legend

cowplot::plot_grid(plotlist = p_list, nrow = 2)

ggsave(
  "viz/Fig4A.multiomics_cor.240217.pdf",
  width = 15, height = 10
)
