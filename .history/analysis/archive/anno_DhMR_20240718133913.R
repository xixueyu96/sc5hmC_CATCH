##---------------------flux ---------------------------
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(ggsankey)
library(data.table)

plot_dt <- readRDS("data/processed/DhMR.sc_5k.rds")

sankey_df <- plot_dt %>%
  # filter(category!="stable") %>%
  mutate(pos=paste0(chr, ":", start, "-", end)) %>%
  mutate(category = ifelse(meth_g2 > 0.5*meth_g1, "active_yes", "active_no")) %>%
  dplyr::select(pos, cmp, category) %>%
  tidyr::spread(key = pos, value = category) %>%
  tibble::column_to_rownames("cmp") %>% t

pal <- c("active_yes"="#ee6300","active_no"="#97ca5b", "NA"="grey")

sankey_long <- sankey_df %>%
  as.data.frame() %>%
  dplyr::select(`LZY->2C`,`2C->4C`, `4C->8C`, `8C->Morula`, `Morula->Blast`) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  group_by(`LZY->2C`,`2C->4C`, `4C->8C`, `8C->Morula`, `Morula->Blast`) %>%
  summarise(comb=n())

## convert NA to no_cover
sankey_long %>%
  filter(comb > 10) %>%
  ggplot(
    aes(
      # axis1 = `Sperm->EZY`,
      axis1 = `LZY->2C`,
      axis2 = `2C->4C`,
      axis3 = `4C->8C`,
      axis4 = `8C->Morula`,
      axis5 = `Morula->Blast`,
      y = comb)) +
  geom_flow(aes(fill=after_stat(stratum))) +
  geom_stratum(aes(fill=after_stat(stratum)), color="white",width = 1/4) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5),
    labels = c("LZY->2C", "2C->4C", "4C->8C", "8C->Morula", "Morula->Blast")) +
  # geom_text(stat = "stratum",
  #           aes(label = after_stat(stratum))) +
  theme_alluvial(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "viz/active_DhMR_flux.5k.2024-05-28.pdf",
  width = 8,
  height = 4
)

head(sankey_df)

consist_active_gr <- sankey_df %>%
  as.data.frame() %>% dplyr::select(-`Sperm->EZY`) %>%
  dplyr::filter(if_all(everything(), ~ . == "active_yes")) %>%
  tibble::rownames_to_column("pos") %>%
  tidyr::separate(
    col = "pos", sep = ":|-",
    into = c("chr", "start", "end")) %>%
  dplyr::select(chr, start, end) %>%
  mutate(region = "consistent_active") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>% sort()

## intersection with hotspot
hotspot_gr <- readBed("data/source/1111_hmC_ratio_1X50.merge.1k.v5s4.hotspot.bed")

res <- ChIPpeakAnno::makeVennDiagram(
  Peaks=list(consist_active_gr,hotspot_gr),
  NameOfPeaks=c("consistent_active", "hotspot")
)

library(Vennerable)
venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}
v <- venn_cnt2venn(res$vennCounts)
plot(v)

## intersection with gene_cor

gene_cor_res <- read.table(
  "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
  header = T, stringsAsFactors = F
)
gene_gr <- readRDS("data/handmade/mm10.genebody.rds")

pos_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho > 0) %>%
  pull(gene)
neg_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho < 0) %>%
  pull(gene)

gene_gr$gene_cor <- case_when(
  gene_gr$gene %in% pos_genes ~ "pos",
  gene_gr$gene %in% neg_genes ~ "neg",
  TRUE ~ "no_sig"
)

res <- ChIPpeakAnno::makeVennDiagram(
  Peaks=list(consist_active_gr, gene_gr[gene_gr$gene_cor=="neg"]),
  NameOfPeaks=c("consistent_active", "neg_cor")
)

v <- venn_cnt2venn(res$vennCounts)
plot(v)


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
        as.data.frame() %>%
        filter(meth_g2 > 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
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

##---------------------active DhMR (homer)---------------------

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

##---------------------active DhMR (ChIPanno)---------------------

if(F){
  library(ChIPseeker)
  library(GenomicRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  stage_list <- c("2C->4C", "4C->8C", "8C->Morula", "Morula->Blast")

  active_AnnoList  <- lapply(

    stage_list,

    function(x){

      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

      in_fn <- file.path(
        "data/processed", paste0(
          gsub("->", "to", x),
          ".DhMR_5mC_ratio.WGBS_GSR.rds"
        )
      )

      tmp_df <- readRDS(in_fn)

      sub_gr <- tmp_df %>%
        filter(meth_g2 < 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      anno_res <- annotatePeak(sub_gr, TxDb=txdb, tssRegion=c(-1000, 1000))

      return(anno_res)
    }
  )

  names(active_AnnoList) <- paste0("(active)", stage_list)

  notactive_AnnoList  <- lapply(

    stage_list,

    function(x){

      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

      in_fn <- file.path(
        "data/processed", paste0(
          gsub("->", "to", x),
          ".DhMR_5mC_ratio.WGBS_GSR.rds"
        )
      )

      tmp_df <- readRDS(in_fn)

      sub_gr <- tmp_df %>%
        filter(meth_g2 < 0.5 * meth_g1) %>%
        dplyr::select(chr, start, end) %>%
        makeGRangesFromDataFrame()

      anno_res <- annotatePeak(sub_gr, TxDb=txdb, tssRegion=c(-1000, 1000))

      return(anno_res)
    }
  )

  names(notactive_AnnoList) <- paste0("(not_active)", stage_list)

  dhmr_AnnoList <- c(active_AnnoList, notactive_AnnoList)

  plotAnnoBar(dhmr_AnnoList)

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


