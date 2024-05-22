##---------------------active DhMR (rGREAT)---------------------

if(F){
  library(rGREAT)
  library(GenomicRanges)
  library(dplyr)

  lapply(

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
        filter(meth_g2 > 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      job <- submitGreatJob(
        sub_gr,species="mm10",
        rule = "basalPlusExt",
        request_interval = 50
      )

      # res_great <- plotRegionGeneAssociationGraphs(job)
      tb_great <- getEnrichmentTables(job)
      # View(tb_great$`GO Biological Process`)
      ## save result
      out_fn <- file.path(
        "data/processed",
        paste0(gsub("->", "to", x),".DhMR_GREAT.tsv"
        ))

      write.table(
        tb_great$`GO Biological Process`, out_fn,
        col.names = T, row.names = F, sep = "\t", quote = F
      )
    }
  )
}

##---------------------active DhMR (LOLA)---------------------

if(F){
  library(dplyr)
  library(circlize)
  library(EnrichedHeatmap)
  library(LOLA)

  ## generate control dataset
  chr_df <-  read.chromInfo(species="mm10")$df
  chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
  chr_window <- makeWindows(chr_gr, w = 5e3)

  ## reference dataset
  regionDB <- loadRegionDB(
    "/mnt/f/packages/LOLACore/mm10/",
    collections = "annotatr"
  )

  lapply(

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
        filter(meth_g2 > 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      userSets <- GRangesList(sub_gr)
      names(userSets) <- x

      ## calc enrichment
      locResults <- runLOLA(
        userSets, chr_window, regionDB,
        minOverlap = 1, direction = "enrichment"
      )

      ## save result
      out_fn <- file.path(
        "data/processed",
        paste0(gsub("->", "to", x),".DhMR_LOLA.tsv"
        ))

      write.table(
        locResults, out_fn,
        col.names = T, row.names = F, sep = "\t", quote = F
      )
    }
  )
}

##---------------------active DhMR (clusterProfiler)---------------------

if(F){
  library(dplyr)
  library(ChIPseeker)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  lapply(

    c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),

    function(x){

      in_fn <- file.path(
        "data/processed", paste0(
          gsub("->", "to", x),
          ".DhMR_5mC_ratio.WGBS_GSR.rds"
        )
      )

      tmp_df <- readRDS(in_fn)

      sub_gr <- tmp_df %>%
        filter(meth_g2 > 0.5 * meth_g1) %>%
        select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      ## calc enrichment
      anno_res <- annotatePeak(sub_gr, TxDb=txdb, tssRegion=c(-1000, 1000))

      entrez <- anno_res@anno %>%
        as.data.frame() %>%
        filter(abs(distanceToTSS) <= 1e3) %>%
        pull(geneId) %>% unique()

      ego <- enrichGO(
        gene = entrez,
        keyType = "ENTREZID",
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
      )

      edox2 <- pairwise_termsim(ego)
      p1 <- treeplot(edox2)
      p2 <- treeplot(edox2, cluster.params = list(method = "ward.D2"))
      aplot::plot_list(p1, p2, tag_levels='A')
    }
  )
}

if(F){
  library(simplifyEnrichment)

  lapply(

    c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),

    function(x){

      in_fn <- file.path(
        "data/processed",
        paste0(gsub("->", "to", x),".DhMR_GREAT.tsv"
      ))

      go_res <- fread(in_fn, header = T, stringsAsFactors = F)

      set.seed(888)
      # go_id <- random_GO(500)
      go_id <- go_res %>%
        filter(Binom_Adjp_BH < 0.01) %>%
        top_n(-500, Binom_Adjp_BH) %>% pull(ID)

      mat <- GO_similarity(go_id, ont = "BP")
      cl <- binary_cut(mat)

      plot_fn <- paste0("viz/", gsub("->", "to", x),".DhMR_GREAT_mclust.pdf")

      pdf(plot_fn, width = 8.5, height = 5)
      ht_clusters(mat, cl)
      dev.off()
    }
  )
}


##---------------------active DhMR (LOLA)---------------------

if(F){
  library(dplyr)
  library(circlize)
  library(EnrichedHeatmap)
  library(LOLA)

  ## generate control dataset
  chr_df <-  read.chromInfo(species="mm10")$df
  chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
  chr_window <- makeWindows(chr_gr, w = 5e3)

  ## reference dataset
  regionDB <- loadRegionDB(
    "/mnt/f/packages/LOLACore/mm10/",
    collections = "annotatr"
  )

  lapply(

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
        filter(meth_g2 > 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      userSets <- GRangesList(sub_gr)
      names(userSets) <- x

      ## calc enrichment
      locResults <- runLOLA(
        userSets, chr_window, regionDB,
        minOverlap = 1, direction = "enrichment"
      )

      ## save result
      out_fn <- file.path(
        "data/processed",
        paste0(gsub("->", "to", x),".DhMR_LOLA.tsv"
        ))

      write.table(
        locResults, out_fn,
        col.names = T, row.names = F, sep = "\t", quote = F
      )
    }
  )
}

if(F){
  library(dplyr)
  library(rtracklayer)

  lapply(

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

      ## test
      sub_gr <- tmp_df %>%
        filter(meth_g2 > 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      out_fn <- file.path(
        "data/processed",
        paste0(gsub("->", "to", x),".homer.bed"
        ))

      write.table(
        sub_gr %>% as.data.frame(),
        file = out_fn, col.names = F, row.names = F, sep = "\t", quote = F
      )

      ## control
      sub_gr <- tmp_df %>%
        filter(!(meth_g2 > 0.5 * meth_g1)) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      out_fn <- file.path(
        "data/processed",
        paste0(gsub("->", "to", x),".homer_control.bed"
        ))

      write.table(
        sub_gr %>% as.data.frame(),
        file = out_fn, col.names = F, row.names = F, sep = "\t", quote = F
      )


      # rtracklayer::export(sub_gr, out_fn, format = "bed")

    }
  )
}

##---------------------active DhMR (clusterProfiler)---------------------

if(F){
  library(dplyr)
  library(ChIPseeker)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  lapply(

    c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"),

    function(x){

      in_fn <- file.path(
        "data/processed", paste0(
          gsub("->", "to", x),
          ".DhMR_5mC_ratio.WGBS_GSR.rds"
        )
      )

      tmp_df <- readRDS(in_fn)

      sub_gr <- tmp_df %>%
        filter(meth_g2 > 0.5 * meth_g1) %>%
        select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      ## calc enrichment
      anno_res <- annotatePeak(sub_gr, TxDb=txdb, tssRegion=c(-1000, 1000))

      entrez <- anno_res@anno %>%
        as.data.frame() %>%
        filter(abs(distanceToTSS) <= 1e3) %>%
        pull(geneId) %>% unique()

      ego <- enrichGO(
        gene = entrez,
        keyType = "ENTREZID",
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
      )

      edox2 <- pairwise_termsim(ego)
      p1 <- treeplot(edox2)
      p2 <- treeplot(edox2, cluster.params = list(method = "ward.D2"))
      aplot::plot_list(p1, p2, tag_levels='A')
    }
  )
}
