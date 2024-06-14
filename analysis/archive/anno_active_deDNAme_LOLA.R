if(F){
  library(dplyr)
  library(circlize)
  library(EnrichedHeatmap)
  library(LOLA)
  library(ggplot2)
  library(magrittr)

  ## reference dataset
  regionDB <- loadRegionDB(
    "/mnt/f/packages/LOLACore/mm10/",
    collections = "annotatr"
  )

  anno_active_DhMR_LOLA <- function(control=c("global", "detected"), save_file=F){
    # control <- "detected"
    # save_file <- F

    anno_res <- lapply(

      c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),

      function(x){

        # stage_pair <- "2C->4C"
        in_fn <- file.path(
          "data/processed", paste0(
            gsub("->", "to", x),
            ".DhMR_5mC_ratio.WGBS_GSR.rds"
          )
        )

        tmp_df <- readRDS(in_fn)

        sub_gr <- tmp_df %>%
          as.data.frame() %>%
          mutate(group = ifelse(meth_g2 > 0.5 * meth_g1, "active", "not_active")) %>%
          # filter(meth_g2 > 0.5 * meth_g1) %>%
          dplyr::select(chr, start, end, group) %>%
          makeGRangesFromDataFrame(keep.extra.columns = T)

        userSets <- GRangesList(
          active = sub_gr[sub_gr$group=="active"],
          not_active = sub_gr[sub_gr$group=="not_active"]
        )
        names(userSets) <- c("active", "not_active")

        ## calc enrichment

        if(control=="global"){

          ## generate control dataset
          chr_df <-  read.chromInfo(species="mm10")$df
          chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
          chr_window <- makeWindows(chr_gr, w = 5e3)

          locResults <- runLOLA(
            userSets, chr_window, regionDB,
            minOverlap = 1, direction = "enrichment") %>%
            mutate(stage_pair = x)
        }
        if(control=="detected"){
          locResults <- runLOLA(
            userSets, sub_gr, regionDB,
            minOverlap = 1, direction = "enrichment")  %>%
            mutate(stage_pair = x)
        }
        # locResults <- runLOLA(
        #   userSets, sub_gr, regionDB,
        #   minOverlap = 1, direction = "enrichment"
        # )

        group_order <- rev(c("genes", "cpg", "enhancer", "encode3Ren", "embryo"))

        locResults %<>%
          # filter(userSet=="active") %>%
          tidyr::separate(
            col = "description",
            sep = "_",
            into = c("annot.group", "annot.type")) %>%
          filter(oddsRatio !=0) %>%
          arrange(factor(annot.group, levels = group_order), oddsRatio)

        locResults$annot.type <- factor(
          locResults$annot.type,
          levels = unique(locResults$annot.type)
        )

        pal <- wesanderson::wes_palette("Zissou1")

        group_order <- rev(c("genes", "cpg", "enhancer", "encode3Ren", "embryo"))

        p <- locResults %>%
          # filter(userSet=="active") %>%
          filter(annot.type!="fantom") %>%
          filter(annot.type!="H3K9me3") %>%
          ggplot(aes(
            x = annot.type,
            # y = log(oddsRatio),
            y = oddsRatio,
            fill = annot.group)) +
          geom_hline(
            yintercept = 1,
            linetype="longdash",
            color="black") +
          geom_bar(
            stat = "identity",
            position = "dodge",
            width = .8,
            alpha = .7)+
          coord_flip() +
          facet_wrap(.~userSet, scales = "free_y") +
          scale_fill_manual(values = c(pal[c(1,3,5)],"#a288ec")) +
          # geom_hline(
          #   yintercept = c(0.5, 1, 1.5),
          #   color="white", size=2) +
          # ylim(c(-2.5, 2.5)) +
          theme_classic()
        print(p)

        if(save_file){
          plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]

          plot_fn <- file.path(
            "viz",
            paste0(gsub("->", "to", x), ".", plot_time,".DhMR_LOLA_annotatr.control_5hmC.pdf"
            ))

          ggsave(
            plot = p, filename = plot_fn, width = 7.5, height = 5.92
          )

          ## save result
          out_fn <- file.path(
            "data/processed",
            paste0(gsub("->", "to", x),".", plot_time,".DhMR_LOLA_annotatr.control_5hmC.tsv"
            ))

          write.table(
            locResults, out_fn,
            col.names = T, row.names = F, sep = "\t", quote = F
          )
        }

        return(locResults)
      }
    )

    anno_res_mrg <- do.call("rbind", anno_res)
    anno_res_mrg %<>%
      filter(oddsRatio > 0) %>%
      filter(annot.group %in% c("cpg", "genes", "encode3Ren")) %>%
      dplyr::select(userSet, annot.group, annot.type, stage_pair, oddsRatio) %>%
      group_by(userSet, annot.group, annot.type) %>%
      reframe(min_odds = min(oddsRatio),
              mean_odds = mean(oddsRatio),
              max_odds = max(oddsRatio)) %>%
      arrange(factor(annot.group, levels = group_order))

    anno_res_mrg$annot.type <- factor(
      anno_res_mrg$annot.type,
      levels = unique(anno_res_mrg$annot.type) %>% rev
    )

    anno_res_mrg$userSet <- factor(
      anno_res_mrg$userSet,
      levels = c("not_active", "active")
    )

    active_pal <- c("not_active"="#9acb60", "active"="#e66200")

    p <- anno_res_mrg %>%
      ggplot(aes(x=annot.type))+
      geom_point(
        aes(y = mean_odds, color = userSet),
        position = position_dodge(width = 2 / 3)) +
      geom_errorbar(
        aes(ymin = min_odds, ymax = max_odds, color = userSet),
        position = position_dodge(width = 2 / 3), width = 0) +
      scale_color_manual(
        values = active_pal,
        name = paste0("ctrl:", control)) +
      # scale_color_discrete(name = paste0("control:", control)) +
      geom_hline(yintercept = 1, linetype="longdash") +
      geom_vline(xintercept = c(6.5, 10.5), linetype="longdash", color="#a9aaad") +
      labs(y= "Odds Ratio") +
      theme_bw(base_size = 15) +
      theme(panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1/3)
    # print(p)
    return(p)
  }

  p1 <- anno_active_DhMR_LOLA(control = "global")
  p2 <- anno_active_DhMR_LOLA(control = "detected")

  cowplot::plot_grid(p1, p2, ncol = 1)

  ggsave(
    "viz/Active_DhMR.anno_LOLA.240614.pdf",
    width = 8.54, height = 5.71, units = "in"
  )
}
