library(LOLA)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

## load input
control_gr <- readRDS("data/processed/hotspot.cpg_density_matched.rds")
hotspot_gr <- readBed("data/source/1111_hmC_ratio_1X50.merge.1k.v5s4.hotspot.bed")
merged_gr <- readRDS("data/processed/hotspot.5hmC.1X50.240218.rds")

## cpgs, genes, enhancer, histone
library(annotatr)
library(rtracklayer)

export_gr <- function(x, fn){
  x <- x %>%
    as.data.frame() %>%
    dplyr::select(
      seqnames, start, end,
      type, strand
    )
  write.table(
    x, fn,
    col.names = F, row.names = F,
    sep = "\t", quote = F
  )
}

mm10_cpg_gr <- build_annotations(genome = 'mm10', annotations = "mm10_cpgs")
mm10_gene_gr <- build_annotations(genome = 'mm10', annotations = "mm10_basicgenes")
mm10_enhancer_gr <- build_annotations(genome = 'mm10', annotations = "mm10_enhancers_fantom")
mm10_intergene_gr <- build_annotations(genome = 'mm10', annotations = "mm10_genes_intergenic")

# R1 mESC line
ab_gr <- import("data/public/4DNFI9685Y2G.bw")
ab_gr <- ab_gr %>% .[!is.na(.$score)]
ab_gr$type <- "A"
ab_gr$type[ab_gr$score < 0] <- "B"
ab_gr <- sort(ab_gr)

export_gr(mm10_cpg_gr, fn = "data/public/mm10.CpGs.bed")
export_gr(mm10_gene_gr, fn = "data/public/mm10.gene.bed")
export_gr(mm10_enhancer_gr, fn = "data/public/mm10.enhancer.bed")
export_gr(mm10_intergene_gr, fn = "data/public/mm10.intergenic.bed")
export_gr(ab_gr, fn = "data/public/mm10.AB_compart.bed")


run_lola <- function(x, control_region=c("matched", "genomic_tile")){

  regionDB <- loadRegionDB(
    "/mnt/f/packages/LOLACore/mm10/",
    collections = x
  )

  userSets <- GRangesList(hotspot_gr)
  names(userSets) <- "hotspot"

  if(control_region=="matched"){

    locResults <- runLOLA(
      userSets,
      # merged_gr,
      c(hotspot_gr, control_gr),
      regionDB,
      minOverlap = 1,
      direction = "enrichment"
    )

  }else if(control_region=="genomic_tile"){
    require(EnrichedHeatmap)

    chr_df <-  circlize::read.chromInfo(species="mm10")$df
    chr_gr <- GRanges(
      seqnames = chr_df[, 1],
      ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3])
    )
    chr_window <- makeWindows(chr_gr, w = unique(width(hotspot_gr)))

    locResults <- runLOLA(
      userSets,
      # merged_gr,
      c(hotspot_gr, chr_window),
      regionDB,
      minOverlap = 1,
      direction = "enrichment"
    )
  }

  plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
  fn_1 <- paste0("data/processed/hotspot.lola_", x, ".", plot_time, ".rds")
  saveRDS(locResults, fn_1)

  p <- locResults %>%
    mutate(label=ifelse(qValue < 0.05, description, "")) %>%
    ggplot(aes(x=oddsRatio, y=-log10(qValue)))+
    geom_point(size=2, shape=21, fill="#f79c56")+
    ggrepel::geom_text_repel(
      aes(label=label),box.padding = 0.5, max.overlaps =5)+
    theme_classic(base_size = 20) +
    theme(aspect.ratio = 1)

  print(p)

  fn_2 <- paste0("viz/hotspot.lola_", x, ".", plot_time, ".pdf")

  # ggsave(
  #   plot = p,
  #   filename = fn_2,
  #   width = 10,
  #   height = 10,
  #   units = "cm"
  # )
}

## repeat
run_lola("rmsk_gene", "matched")
run_lola("rmsk_LINE", "matched")
run_lola("rmsk_LTR", "matched")


## genomic elements
run_lola("annotatr", "matched")

locResults <- readRDS("data/processed/hotspot.lola_annotatr.2024-02-20.rds")

pal <- wesanderson::wes_palette("Zissou1")

group_order <- rev(c("genes", "cpg", "enhancer", "encode3Ren", "embryo"))

locResults %<>%
  tidyr::separate(
    col = "description",
    sep = "_",
    into = c("annot.group", "annot.type")) %>%
  filter(oddsRatio !=0) %>%
  arrange(factor(annot.group, levels = group_order), oddsRatio)

locResults$annot.type <- factor(
  locResults$annot.type,
  levels = locResults$annot.type
)

locResults %>%
  filter(annot.type!="fantom") %>%
  filter(annot.type!="H3K9me3") %>%
  ggplot(aes(
    x = annot.type,
    # y = log(oddsRatio),
    y = oddsRatio,
    fill = annot.group)) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    width = .8,
    alpha = .7)+
  coord_flip() +
  scale_fill_manual(values = c(pal[c(1,3,5)],"#a288ec")) +
  # geom_hline(
  #   yintercept = c(0.5, 1, 1.5),
  #   color="white", size=2) +
  # ylim(c(-2.5, 2.5)) +
  geom_hline(
    yintercept = 1,
    linetype="longdash",
    color="black") +
  theme_classic()

ggsave(
  "viz/hotspot.lola_annotatr.2024-02-20.pdf",
  width = 5, height = 6
)

## H3K9me3
h3k9me3_gr <- LOLA::readBed("data/processed/H3K9me3_GSR.broadPeak.stable.bed")

hotspot_hit <- findOverlaps(hotspot_gr, h3k9me3_gr)
control_hit <- findOverlaps(control_gr, h3k9me3_gr)

chi_input <- matrix(c(
  length(hotspot_hit@from),
  length(control_hit@from),
  length(hotspot_gr) - length(hotspot_hit@from),
  length(control_gr) - length(control_hit@from)
), nrow = 2)
rownames(chi_input) <- c("hotspot_hit", "control")
colnames(chi_input) <- c("hit", "not_hit")

chisq.test(chi_input)

chi_input[1,1]*chi_input[2,2]/(chi_input[1,2]*chi_input[2,1])




