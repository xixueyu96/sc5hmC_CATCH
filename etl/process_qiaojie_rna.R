library(data.table)
library(org.Mm.eg.db)
library(GenomicFeatures)

##### prepare annotation #####

if(F){

  ## download gtf from GENCODE
  # download.file(
  #   "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/gencode.vM11.annotation.gtf.gz",
  #   "data/public/gencode.vM11.annotation.gtf.gz",
  #   method = "wget"
  # )

  txdb <- makeTxDbFromGFF("data/public/gencode.vM11.annotation.gtf",format="gtf")
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
  gene_length <- data.frame(exonic.gene.sizes)
  gene_length$ensembl <- sapply(strsplit(rownames(gene_length),split='[.]'),function(x) x[1])

  gene_length$symbol <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    gene_length$ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = 'first'
  )

  saveRDS(gene_length, file = 'data/public/gencode.vM11.gene_length.Rds')

}

##### data cleaning #####
if (F) {
  rna.mtx <-
    fread(
      'data/public/GSE136714_raw.counts.txt'
    ) %>% as.data.frame()

  rownames(rna.mtx) <- rna.mtx[,1]
  rna.mtx <- rna.mtx[, -1]

  gmap <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    rownames(rna.mtx),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = 'first'
  )

  index <- is.na(gmap)
  gmap[index] <- names(gmap)[index]

  # aggregate by gene ID
  rep.id <- names(gmap)[gmap %in% names(which(table(gmap) > 1))]

  rna.mtx.unique <- rna.mtx[!rownames(rna.mtx) %in% rep.id, ] %>% as.matrix()
  rownames(rna.mtx.unique) <- gmap[rownames(rna.mtx.unique)]

  agg.mtx <- aggregate(rna.mtx[rep.id, ], list(gene=gmap[rep.id]),sum)
  rownames(agg.mtx) <- agg.mtx$gene
  agg.mtx <- agg.mtx[,-1] %>% as.matrix

  rna.mtx <- rbind(rna.mtx.unique,agg.mtx)

  # aggregate split cells
  agg.mtx <- aggregate(t(rna.mtx),list(cell.id=colnames(rna.mtx) %>% gsub('_split.','',.)),mean)
  rownames(agg.mtx) <- agg.mtx$cell.id
  agg.mtx <- agg.mtx[,-1] %>% as.matrix
  rna.mtx <- t(agg.mtx)

  saveRDS(rna.mtx, 'data/processed/GSE136714.aggregated_gene_matrix.Rds')

}
