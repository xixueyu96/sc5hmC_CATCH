
plot_DHMR.frac.v1 <- function(plot_dt, save_file=F, level=input_level){

  require(data.table)
  require(dplyr)
  require(scales)
  require(ggplot2)

  plot_dt$cmp <- factor(plot_dt$cmp, levels = level)

  p <- plot_dt %>%
    filter(cmp %in% level) %>%
    mutate(type=case_when(
      category=="stable" ~ "stable",
      category!="stable" ~ "changing"
    )) %>%
    group_by(cmp, type) %>%
    summarise(tile=n()) %>%
    ggplot(aes(x=cmp, y=tile,fill=type)) +
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(y="Fraction of tiles",
         title = paste0("tile_width:", unique(plot_dt$tile_width)))+
    scale_fill_manual(values = type_col_list)+
    theme_classic(base_size = 20)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if(save_file){
    fileName <- paste("../plot/DMHR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.frac_stable.pdf",sep = ".")
    ggsave(p, filename = fileName)
  } else print(p)
}
plot_DHMR.frac.v2 <- function(plot_dt, save_file=F, level=input_level){

  require(data.table)
  require(dplyr)
  require(scales)
  require(ggplot2)

  plot_dt$cmp <- factor(plot_dt$cmp, levels = level)

  p <- plot_dt %>%
    filter(category!="stable") %>%
    group_by(cmp, category) %>%
    summarise(tile=n()) %>%
    ggplot(aes(x=cmp, y=tile,fill=category)) +
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(y="Fraction of tiles",
         title = paste0("tile_width:", unique(plot_dt$tile_width)))+
    scale_fill_manual(values = category_col_list)+
    theme_classic(base_size = 20)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if(save_file){
    fileName <- paste("../plot/DMHR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.frac_changing.pdf",sep = ".")
    ggsave(p, filename = fileName)
  } else print(p)
}
plot_DHMR.frac.stable <- function(plot_dt, save_file=F, level=input_level){

  require(data.table)
  require(dplyr)
  require(scales)
  require(ggplot2)

  plot_dt$cmp <- factor(plot_dt$cmp, levels = level)

  p <- plot_dt %>%
    filter(cmp %in% level) %>%
    mutate(type=case_when(
      category=="stable" ~ "stable",
      category!="stable" ~ "changing"
    )) %>%
    group_by(cmp, type) %>%
    summarise(tile=n()) %>%
    ggplot(aes(x=cmp, y=tile,fill=type)) +
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(title = "Fraction of tiles")+
    scale_fill_manual(values = type_col_list)+
    coord_flip()+
    scale_y_continuous(
      position = "right")+
    theme_classic(base_size = 20)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size=10),
          legend.direction = "vertical",
          legend.position = "bottom")

  if(save_file){
    fileName <- paste("../plot/DMHR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.frac_stable.pdf",sep = ".")
    ggsave(p, filename = fileName)
  } else return(p)
}
plot_DHMR.frac.changing <- function(plot_dt, save_file=F, level=input_level){

  require(data.table)
  require(dplyr)
  require(scales)
  require(ggplot2)

  plot_df <- plot_dt %>% as.data.frame()
  plot_df$cmp <- factor(plot_df$cmp, levels = level)
  plot_df$category <- factor(plot_df$category, levels = c("stable", "increasing", "decreasing"))

  p <- plot_df %>%
    filter(cmp %in% level) %>%
    dplyr::group_by(cmp, category, .drop=FALSE) %>%
    dplyr::summarise(tile=n()) %>%
    filter(category!="stable") %>%
    ggplot(aes(x=cmp, y=tile,fill=category)) +
    geom_bar(position="fill", stat="identity", width = .8)+
    labs(title = "Fraction of tiles")+
    scale_fill_manual(values = category_col_list)+
    coord_flip()+
    scale_y_continuous(
      position = "right")+
    theme_classic(base_size = 20)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size=10),
          legend.direction = "vertical",
          legend.position = "bottom")

  if(save_file){
    fileName <- paste("../plot/DMHR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.frac_changing.pdf",sep = ".")
    ggsave(p, filename = fileName)
  } else return(p)
}
plot_DHMR.abs.stable <- function(plot_dt, save_file=F, level=input_level){

  require(data.table)
  require(dplyr)
  require(scales)
  require(ggplot2)

  plot_dt$cmp <- factor(plot_dt$cmp, levels = level)

  p <- plot_dt %>%
    filter(cmp %in% level) %>%
    mutate(type=case_when(
      category=="stable" ~ "stable",
      category!="stable" ~ "changing"
    )) %>%
    group_by(cmp, type) %>%
    summarise(tile=n()) %>%
    ggplot(aes(x=cmp, y=tile,fill=type)) +
    geom_bar(stat="identity", width = .8)+
    labs(title = "Absolute number of tiles")+
    scale_fill_manual(values = type_col_list)+
    coord_flip()+
    scale_y_continuous(
      position = "right",labels = label_scientific())+
    theme_classic(base_size = 20)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size=10),
          legend.direction = "vertical",
          legend.position = "bottom")

  if(save_file){
    fileName <- paste("../plot/DMHR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.frac_changing.pdf",sep = ".")
    ggsave(p, filename = fileName)
  } else print(p)
}
plot_DHMR.abs.changing <- function(plot_dt, save_file=F, level=input_level){

  require(data.table)
  require(dplyr)
  require(scales)
  require(ggplot2)

  plot_dt$cmp <- factor(plot_dt$cmp, levels = level)
  plot_dt$category <- factor(plot_dt$category, levels = c("stable", "increasing", "decreasing"))

  p <- plot_dt %>%
    filter(cmp %in% level) %>%
    group_by(cmp, category, .drop = FALSE) %>%
    summarise(tile=n()) %>%
    filter(category!="stable") %>%
    ggplot(aes(x=cmp, y=tile,fill=category)) +
    geom_bar(stat="identity", width = .8)+
    # labs(x="Absolute number of tiles",
    #      title = paste0("tile_width:", unique(plot_dt$tile_width)))+
    labs(title = "Absolute number of tiles")+
    scale_fill_manual(values = category_col_list)+
    coord_flip()+
    scale_y_continuous(
      position = "right",labels = label_scientific())+
    theme_classic(base_size = 20)+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          plot.title = element_text(size=10),
          legend.direction = "vertical",
          legend.position = "bottom")

  if(save_file){
    fileName <- paste("../plot/DMHR", unique(plot_dt$tile_width),
                      "qval0.05.minLFC_1.frac_changing.pdf",sep = ".")
    ggsave(p, filename = fileName)
  } else print(p)
}
