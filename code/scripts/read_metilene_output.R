library(data.table)
library(dplyr)
library(scales)

cmp_list <- c("Sperm","LZY","2C","4C","8C","Morula")
names(cmp_list) <- c("EZY","2C","4C","8C","Morula", "Blast")
levels <- c(
  "Sperm->EZY","LZY->2C",
  "2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"
)

merge_DHMR_lfc <- function(compare_list,tile_width, min_lfc=1){
  
  todo_index <- 1:length(compare_list)
  
  worker_fun <- function(x){
    
    pre_stage <- compare_list[x]
    pro_stage <- names(compare_list[x])
    
    ## input
    in_dir <- "F:/Gitee/sc5hmC/data/processed/"
    fileName <- paste0("Stage_",pre_stage,"_to_",pro_stage,".",tile_width,".DMR.metilene.out.gz")
    fileName <- file.path(in_dir, fileName)
    
    if(file.info(fileName)$size==0) return(NULL)
    else{
      dt <- fread(fileName, stringsAsFactors = F)
      dt <- dt[,c(1:5,9:10)]
      colnames(dt) <- c("chr", "start", "end", "qval", "meth_diff", "meth_g1", "meth_g2")
      
      dt_summary <- dt %>% filter(meth_g1 >=0.001 & meth_g2 >=0.001) %>%
        mutate(lfc_meth=log2(meth_g2/meth_g1),
               cmp=paste0(pre_stage, "->", pro_stage)) %>%
        mutate(category=case_when(
          qval < 0.05 & lfc_meth >= min_lfc ~ "increasing",
          qval >=0.05 | (qval < 0.05 & abs(lfc_meth) < min_lfc ) ~"stable",
          qval < 0.05 & lfc_meth <= -min_lfc ~"decreasing"
        ))
      
      return(dt_summary)
    }
    # stopifnot(file.info(fileName)$size!=0)
  }
  
  res_list <- lapply(todo_index, worker_fun)
  # res_list <- res_list[-which(sapply(res_list, is.null))]
  res <- do.call("rbind", res_list)
  res$tile_width <- tile_width
  
  return(res)
}

plot_dt <- merge_DHMR_lfc(cmp_list, "5k", min_lfc = 1)
