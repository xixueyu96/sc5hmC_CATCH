library(data.table)
library(dplyr)

L2C_rna <- fread("data/public/GSE135457_total_RNA_FPKM_all.txt.gz")
E2C_rna <- fread("data/public/GSM4476421_2cell_early_total_RNA_FPKM.txt.gz")
mii_rna <- fread("data/public/GSM4476422_MII_total_RNA_FPKM.txt.gz")

mii_gene <- mii_rna %>%
  filter(MII_rep1> 5) %>%
  pull(gene)

zga_major <- intersect(
  L2C_rna %>%
    filter(`2cell_late_con` >5) %>%
    pull(gene),
  mii_rna %>%
    filter(MII_rep1< 5) %>%
    pull(gene)
)
zga_1C_minor <- intersect(
  L2C_rna %>%
    filter(`1cell` >2) %>%
    pull(gene),
  mii_rna %>%
    filter(MII_rep1< 0.5) %>%
    pull(gene)
)
zga_2C_minor <- intersect(
  E2C_rna %>%
    filter(`2cell_early` >2) %>%
    pull(gene),
  mii_rna %>%
    filter(MII_rep1< 0.5) %>%
    pull(gene)
)

library(ggvenn)

ggvenn::ggvenn(
  list(
    "major ZGA" = zga_major,
    "1C minor ZGA" = zga_1C_minor,
    "2C minor ZGA"=zga_2C_minor
    )
)

gene_cor_res <- read.table(
  "data/processed/Gene_wise_cor.5hmC_5mC.20240217.tsv",
  header = T, stringsAsFactors = F
)

pos_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho > 0) %>%
  arrange(rho)

neg_genes <- gene_cor_res %>%
  filter(pval < 0.05 & rho < 0) %>%
  arrange(-rho)

ggvenn(
  list(
    "ZGA_major" = zga_major,
    "pos_gene" = pos_genes$gene,
    "neg_gene" = neg_genes$gene
  )
)

minor_major_common <- intersect(
  zga_major,c(zga_1C_minor, zga_2C_minor)
)

zga_major <- setdiff(zga_major, minor_major_common)
zga_minor <- setdiff(c(zga_1C_minor, zga_2C_minor), minor_major_common)

lt <- list(
  "ZGA_major" = zga_major,
  "ZGA_minor" = zga_minor,
  # "ZGA_1C_minor" = zga_1C_minor,
  # "ZGA_2C_minor" = zga_2C_minor,
  "pos_gene" = pos_genes$gene,
  "neg_gene" = neg_genes$gene
)

m <- make_comb_mat(lt)
UpSet(
  m[comb_size(m) >= 4],
  top_annotation = upset_top_annotation(m[comb_size(m) >= 4], add_numbers = TRUE),
  set_order = c("pos_gene", "neg_gene", "ZGA_major", "ZGA_minor"),
  comb_order = order(comb_size(m)[comb_size(m) >= 4], decreasing = T),
)



