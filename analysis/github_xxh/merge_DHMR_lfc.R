
merge_DHMR_lfc <- function(indir,compare_list,tile_width, min_lfc=1){

  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(scales)
  require(egg)
  require(RColorBrewer)

  todo_index <- 1:length(compare_list)

  worker_fun <- function(x){

    pre_stage <- compare_list[x]
    pro_stage <- names(compare_list[x])

    ## input
    # fileName <- paste0("Stage_",pre_stage,"_to_",pro_stage,".",tile_width,"_qval.0.05.out")
    fileName <- paste0(indir,"Stage_",pre_stage,"_to_",pro_stage,".",tile_width,".DMR.metilene.out")

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
