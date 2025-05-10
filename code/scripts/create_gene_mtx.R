library(data.table)
library(dplyr)
library(Biobase)

##--------------------5hmC--------------------
hmc_mtx <- lapply(
  list.files("data/source/hmC_genebody/", full.names = T),
  function(x){
    df <- read.table(x, header = F, stringsAsFactors = F)
    sp <- strsplit(basename(x), "[.]")[[1]][1]
    # sp <- strsplit(sp, "_")[[1]][2]
    colnames(df) <- c("chr", "start", "end", "gene", sp)
    df <- df %>% .[,c("gene", sp)]
    print(sapply(df, class))
    return(df)
  }
)
hmc_mtx <- hmc_mtx %>% purrr::reduce(function(x, y)
  base::merge(x, y, by = c("gene"))) %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.))

colnames(hmc_mtx) <-
  c("gene",
    "Two-cell",
    "Four-cell",
    "Eight-cell",
    "Blast",
    "Morula",
    "Oocyte",
    "Sperm",
    "Zygote")

head(hmc_mtx)

write.table(
  hmc_mtx,
  file = "data/processed/hmC.genebody_mtx.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F
)

##--------------------5mC--------------------
## methylC-seq data of Liu Jiang
mc_mtx_lj <-
  read.table(
    "data/source/mC_gene_ratio.tsv",
    header = T,
    stringsAsFactors = F
  )
mc_mtx_lj <- mc_mtx_lj %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.)) %>%
  dplyr::select(gene, oocyte, sperm)

mc_mtx_gsr <- lapply(
  list.files("data/source/mC_genebody/", full.names = T),
  function(x){
    df <- read.table(x, header = F, stringsAsFactors = F)
    sp <- strsplit(basename(x), "[.]")[[1]][1]
    sp <- strsplit(sp, "_")[[1]][2]
    colnames(df) <- c("chr", "start", "end", "gene", sp)
    df <- df %>% .[,c("gene", sp)]
    print(sapply(df, class))
    return(df)
  }
)

mc_mtx_gsr <-
  mc_mtx_gsr %>% purrr::reduce(function(x, y)
    base::merge(x, y, by = c("gene"))) %>%
  filter(isUnique(gene)) %>%
  filter(complete.cases(.))

mc_mtx <- inner_join(
  mc_mtx_lj, mc_mtx_gsr,by = "gene")

colnames(mc_mtx) <-
  c("gene",
    "Oocyte",
    "Sperm",
    "Zygote",
    "Two-cell",
    "Four-cell",
    "Eight-cell",
    "Morula",
    "Blast")

head(mc_mtx)

write.table(
  mc_mtx,
  file = "data/processed/mC.genebody_mtx.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F
)

##--------------------RNA--------------------
library(Seurat)
library(egg)

mouse_preim <- readRDS("data/processed/Mouse_preimplan.seurat.rds")
rna_mtx <- (mouse_preim@assays$RNA@data %>% as.matrix() %>% exp)-1

rna_mtx <- rna_mtx %>%
  reshape2::melt() %>%
  left_join(
    mouse_preim@meta.data %>%
      select(c("Cell_id", "Stage")),
    by = c("Var2" = "Cell_id")) %>%
  mutate(
    Stage = case_when(
      Stage %in% c("E2C", "L2C") ~ "2C",
      Stage %in% c("ICM", "TE") ~ "Blast",
      TRUE ~ Stage )) %>%
  dplyr::group_by(Stage, Var1) %>%
  # dplyr::summarise(mean_expr=mean(value, na.rm=T)) %>%
  dplyr::summarise(mean_expr=log1p(mean(value, na.rm=T))) %>%
  tidyr::spread(key = Stage, value = mean_expr)

rna_mtx <- rna_mtx %>%
  dplyr::select(
    Var1, `M2`, `ZY`,
    `2C`, `4C`, `8C`,
    `Morula`, `Blast`
  )
colnames(rna_mtx) <-
  c("gene",
    "Oocyte",
    "Zygote",
    "Two-cell",
    "Four-cell",
    "Eight-cell",
    "Morula",
    "Blast")

head(rna_mtx)

write.table(
  rna_mtx,
  file = "data/processed/RNA.genebody_mtx.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F
)

