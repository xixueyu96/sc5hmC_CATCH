
source("etl/merge_DHMR_lfc.R")
source("etl/plot_DHMR.frac.R")
source("etl/plot_RDL_RhML.R")

type_col_list <- c("stable"="#a9caec","changing"="#4a7bb7")
category_col_list <- c("increasing"="#f16c78","decreasing"="#74bbe7")

cmp_list <- c("Sperm","Oocyte","EZY","LZY","2C","4C","8C","Morula")
names(cmp_list) <- c("EZY","EZY","LZY","2C","4C","8C","Morula", "Blast")
levels <- c(
  "Sperm->EZY","Oocyte->EZY", "EZY->LZY","LZY->2C",
  "2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"
)

plot_dt <- merge_DHMR_lfc("data/source/metilene_out/",cmp_list, "5k", min_lfc = 1)
p2 <- plot_DHMR.frac.stable(plot_dt, save_file=F, level=rev(levels))
p2 <- p2+theme(axis.text.y = element_blank())
p1 <- plot_DHMR.abs.stable(plot_dt, save_file=F, level=rev(levels))
p4 <- plot_DHMR.frac.changing(plot_dt, save_file=F, level=rev(levels))
p4 <- p4+theme(axis.text.y = element_blank())
p3 <- plot_DHMR.abs.changing(plot_dt, save_file=F, level=rev(levels))
p3 <- p3+theme(axis.text.y = element_blank())
p5 <- plot_RDL.v2(plot_dt, save_file=F, level=rev(levels))
p5 <- p5+theme(axis.text.y = element_blank())
p6 <- plot_RhML.v2(plot_dt, save_file=F, level=rev(levels))
p6 <- p6+theme(axis.text.y = element_blank())
egg::ggarrange(p1,p2,p3,p4,p5,p6,nrow = 1)


