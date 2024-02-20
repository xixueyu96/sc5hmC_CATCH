
output_DHMR <- function(plot_dt, tile_width){

  options(scipen = 200) ## scancel scientific notation

  twd_list <- c(1e3, 5e3, 1e4, 5e4)
  names(twd_list) <- c("1k", "5k", "10k", "50k")

  twd <- as.numeric(twd_list[tile_width])

  for (i in names(table(plot_dt$cmp))) {
    for(j in names(table(plot_dt$category))){
      out_df <- plot_dt[plot_dt$cmp==i & plot_dt$category==j,]
      out_df$start <- out_df$end-twd
      out_df <- as.data.frame(out_df)

      # fileName <- paste("../DHMR", gsub("->", "_", i), tile_width, "qval0.05.minLFC_1.bed",sep = ".")
      fileName <- paste("../DHMR", gsub("->", "_", i), tile_width,j, "qval0.05.thresh_0.5.bed",sep = ".")
      write.table(out_df, fileName, col.names = F, row.names = F, sep = "\t", quote = F)
      print(paste0(i, "_", j, " is done!"))
    }

  }

}
