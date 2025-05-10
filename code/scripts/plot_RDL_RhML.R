
plot_RDL <- function(plot_dt, save_file=T, level=input_level){

  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(scales)

  rdl_dt <- plot_dt %>%
    filter(category=="decreasing") %>%
    # mutate(RDL=(meth_g1-meth_g2)/meth_g1) %>%
    mutate(RDL=meth_g2/meth_g1) %>%
    mutate(RDL_bin=cut(RDL, breaks=seq(0,1,0.1), include.lowest=T)) %>%
    group_by(cmp,RDL_bin) %>%
    summarise(RDL_bin_n =n())

  rdl_dt$cmp <- factor(rdl_dt$cmp, levels = level)

  col_list <- colorRampPalette(brewer.pal(n = 5, name =  "YlGnBu"))(length(unique(rdl_dt$RDL_bin)))
  p <- rdl_dt %>%
    ggplot(aes(x=" ", y=RDL_bin_n, fill=RDL_bin))+
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(y="Fraction of binned RDL")+
    scale_fill_manual(values = col_list)+
    coord_polar("y", start=0)+
    facet_grid(.~cmp)+
    theme_void(base_size = 20)+
    theme(plot.title = element_text(hjust = 0.5))
  # theme_classic(base_size = 20)+
  # theme(axis.title.x = element_blank(),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.text.x = element_text(angle = 45, hjust = 1))
  if(save_file){
    fileName <- paste("../plot/DhMR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.RDL_bin_frac.pdf",sep = ".")
    ggsave(p, filename = fileName,width = length(unique(rdl_dt$cmp))*3.75, height = 4)
  } else return(p)
}
plot_RhML <- function(plot_dt, save_file=T, level=input_level){

  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(scales)

  rdl_dt <- plot_dt %>%
    filter(category=="increasing") %>%
    # mutate(RDL=(meth_g2-meth_g1)/meth_g2) %>%
    mutate(RDL=meth_g1/meth_g2) %>%
    mutate(RhML_bin=cut(RDL, breaks=seq(0,2,0.25), include.lowest=T)) %>%
    group_by(cmp,RhML_bin) %>%
    summarise(RDL_bin_n =n())

  rdl_dt$cmp <- factor(rdl_dt$cmp,level)

  col_list <- colorRampPalette(brewer.pal(n = 5, name =  "YlOrRd"))(length(unique(rdl_dt$RhML_bin)))
  p <- rdl_dt %>%
    ggplot(aes(x=" ", y=RDL_bin_n, fill=RhML_bin))+
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(y="Fraction of binned RhML")+
    scale_fill_manual(values = col_list)+
    coord_polar("y", start=0)+
    facet_grid(.~cmp)+
    theme_void(base_size = 20)+
    theme(plot.title = element_text(hjust = 0.5))
  # theme_classic(base_size = 20)+
  # theme(axis.title.x = element_blank(),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.text.x = element_text(angle = 45, hjust = 1))
  if(save_file){
    fileName <- paste("../plot/DhMR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.RhML_bin_frac.pdf",sep = ".")
    ggsave(p, filename = fileName,
           width = length(unique(rdl_dt$cmp))*3.75, height = 4)
  } else return(p)
}
plot_RDL.v2 <- function(plot_dt, save_file=T, level=input_level){

  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(scales)
  require(RColorBrewer)

  rdl_dt <- plot_dt %>%
    # filter(category=="decreasing") %>%
    filter(meth_g1-meth_g2 > 0) %>%
    mutate(RDL=(meth_g1-meth_g2)/meth_g1) %>%
    # mutate(RDL=meth_g2/meth_g1) %>%
    mutate(RDL_bin=cut(RDL, breaks=seq(0,1,0.2), include.lowest=T)) %>%
    group_by(cmp,RDL_bin) %>%
    summarise(RDL_bin_n =n())

  rdl_dt$cmp <- factor(rdl_dt$cmp, levels = level)

  col_list <- colorRampPalette(brewer.pal(n = 3, name =  "YlGnBu"))(length(unique(rdl_dt$RDL_bin)))
  # p <- rdl_dt %>%
  #   ggplot(aes(x=" ", y=RDL_bin_n, fill=RDL_bin))+
  #   geom_bar(position="fill", stat="identity", width = .8)+
  #   labs(y="Fraction of binned RDL")+
  #   scale_fill_manual(values = col_list)+
  #   coord_polar("y", start=0)+
  #   facet_wrap(.~cmp, ncol = 1)+
  #   theme_void(base_size = 10)+
  #   theme(plot.title = element_blank(),
  #         strip.text = element_blank(),
  #         legend.direction = "vertical",
  #         legend.position = "bottom")+
  #   theme(aspect.ratio = 1)
  p <- rdl_dt %>%
    ggplot(aes(x=cmp, y=RDL_bin_n, fill=RDL_bin))+
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(title="Fraction of binned RDL")+
    scale_fill_manual(values = col_list)+
    coord_flip()+
    scale_y_continuous(
      position = "right",labels = label_number())+
    theme_classic(base_size = 20)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size=10),
          legend.direction = "vertical",
          legend.position = "bottom")
  # theme_classic(base_size = 20)+
  # theme(axis.title.x = element_blank(),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.text.x = element_text(angle = 45, hjust = 1))
  if(save_file){
    fileName <- paste("../plot/DhMR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.RDL_bin_frac.pdf",sep = ".")
    ggsave(p, filename = fileName,width = length(unique(rdl_dt$cmp))*3.75, height = 4)
  } else return(p)
}
plot_RhML.v2 <- function(plot_dt, save_file=T, level=input_level){

  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(scales)
  require(RColorBrewer)

  rdl_dt <- plot_dt %>%
    # filter(category=="increasing") %>%
    filter(meth_g2-meth_g1 > 0) %>%
    mutate(RDL=(meth_g2-meth_g1)/meth_g2) %>%
    # mutate(RDL=meth_g1/meth_g2) %>%
    mutate(RhML_bin=cut(RDL, breaks=seq(0,1,0.2), include.lowest=T)) %>%
    group_by(cmp,RhML_bin) %>%
    summarise(RDL_bin_n =n())

  rdl_dt$cmp <- factor(rdl_dt$cmp,level)

  col_list <- colorRampPalette(brewer.pal(n = 3, name =  "YlOrRd"))(length(unique(rdl_dt$RhML_bin)))
  # p <- rdl_dt %>%
  #   ggplot(aes(x=" ", y=RDL_bin_n, fill=RhML_bin))+
  #   geom_bar(position="fill", stat="identity", width = .8)+
  #   labs(y="Fraction of binned RhML")+
  #   scale_fill_manual(values = col_list)+
  #   coord_polar("y", start=0)+
  #   facet_wrap(.~cmp, ncol = 1)+
  #   theme_void(base_size = 10)+
  #   theme(plot.title = element_blank(),
  #         strip.text = element_blank(),
  #         legend.direction = "vertical",
  #         legend.position = "bottom")+
  #   theme(aspect.ratio = 1)
  p <- rdl_dt %>%
    ggplot(aes(x=cmp, y=RDL_bin_n, fill=RhML_bin))+
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(title="Fraction of binned RhML")+
    scale_fill_manual(values = col_list)+
    coord_flip()+
    scale_y_continuous(
      position = "right",labels = label_number())+
    theme_classic(base_size = 20)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size=10),
          legend.direction = "vertical",
          legend.position = "bottom")
  # theme_classic(base_size = 20)+
  # theme(axis.title.x = element_blank(),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.text.x = element_text(angle = 45, hjust = 1))
  if(save_file){
    fileName <- paste("../plot/DhMR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.RhML_bin_frac.pdf",sep = ".")
    ggsave(p, filename = fileName,
           width = length(unique(rdl_dt$cmp))*3.75, height = 4)
  } else return(p)
}
