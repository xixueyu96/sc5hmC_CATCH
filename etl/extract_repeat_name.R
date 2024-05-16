library(data.table)
library(dplyr)

ucsc_rmsk_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz"
desc_dir <- "data/public/mm10.rmsk.txt.gz"

download.file(ucsc_rmsk_url,desc_dir)

total_rmsk <- fread("data/public/mm10.rmsk.txt.gz")

repeat_name <- list(
  L1 = total_rmsk %>% filter(V13=="L1") %>%
    dplyr::select(V11) %>% distinct() %>% pull(V11),
  L2 = total_rmsk %>% filter(V13=="L2") %>%
    dplyr::select(V11) %>% distinct() %>% pull(V11),
  Alu = total_rmsk %>% filter(V13=="Alu") %>%
    dplyr::select(V11) %>% distinct() %>% pull(V11),
  MIR = total_rmsk %>% filter(V13=="MIR") %>%
    dplyr::select(V11) %>% distinct() %>% pull(V11),
  ERVK = total_rmsk %>% filter(V13=="ERVK") %>%
    dplyr::select(V11) %>% distinct() %>% pull(V11),
  ERVL = total_rmsk %>% filter(V13=="ERVL") %>%
    dplyr::select(V11) %>% distinct() %>% pull(V11)
)

saveRDS(repeat_name, "data/processed/mm10.repeat_name.rds")
