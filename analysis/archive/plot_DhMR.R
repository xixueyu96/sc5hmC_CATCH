
source("etl/merge_DHMR_lfc.R")
source("etl/plot_DHMR.frac.R")
source("etl/plot_RDL_RhML.R")

type_col_list <- c("stable"="#a9caec","changing"="#4a7bb7")
category_col_list <- c("increasing"="#f16c78","decreasing"="#74bbe7")

cmp_list <- c("Sperm","LZY","2C","4C","8C","Morula")
names(cmp_list) <- c("EZY","2C","4C","8C","Morula", "Blast")
levels <- c(
  "Sperm->EZY","LZY->2C",
  "2C->4C", "4C->8C", "8C->Morula", "Morula->Blast"
)

plot_dt <- readRDS("data/processed/DhMR.sc_5k.rds")
# plot_dt <- merge_DHMR_lfc("data/source/metilene_out/",cmp_list, "5k", min_lfc = 1)
p2 <- plot_DHMR.frac.stable(plot_dt, save_file=F, level=rev(levels))
p2 <- p2+theme(axis.text.y = element_blank())
p1 <- plot_DHMR.abs.stable(plot_dt, save_file=F, level=rev(levels))
p4 <- plot_DHMR.frac.changing(plot_dt, save_file=F, level=rev(levels))
p4 <- p4+theme(axis.text.y = element_blank())
p3 <- plot_DHMR.abs.changing(plot_dt, save_file=F, level=rev(levels))
p3 <- p3+theme(axis.text.y = element_blank())
p5 <- plot_RDL.v2(plot_dt, save_file=F, level=rev(levels))
p5 <- p5+theme(axis.text.y = element_blank())
p6 <- plot_RhML.v2(plot_dt, save_file=F, level=rev(levels))
p6 <- p6+theme(axis.text.y = element_blank())
egg::ggarrange(p1,p2,p3,p4,p5,p6,nrow = 1)

saveRDS(plot_dt, file = "data/processed/DhMR.sc_5k.rds")

LZto2C <- plot_dt %>% filter(cmp=="LZY->2C") %>%
  mutate(start = end - 5000) %>%
  select(chr, start, end, category)

## plot 5hmC abundance in WT and TETcKO late zygote and 2-cell
options(scipen=999)
write.table(
  LZto2C, file = "data/processed/LZYto2C.DhmR.sc_5k.bed",
  quote = F, sep = "\t", col.names = F, row.names = F
)
options(scipen=2)
# bash calc_5hmC_in_DhmR.sh WT_2C.sample.xls

meta <- lapply(
  list.files("data/documentation", pattern = "sample.xls", full.names = T),
  function(x){
    df <- read.table(x) %>% setNames("Sample") %>%
      mutate(
        Group = gsub(".sample.xls", "", basename(x))
      )
  }
)

meta <- do.call("rbind", meta)

DhmR_5hmC_abundance <- read.table(
    "data/processed/WT_TETcKO.LZY_to_2C.DhmR.sc_5k.abundance.tsv",
    header = F, stringsAsFactors = F) %>%
  setNames(c("Sample", "category", "abundance")) %>%
  inner_join(meta, by = "Sample")

DhmR_5hmC_abundance$category <- factor(
  DhmR_5hmC_abundance$category,
  levels = c("stable", "decreasing", "increasing")
)

DhmR_5hmC_abundance$Group <- factor(
  DhmR_5hmC_abundance$Group,
  levels = c("WT_LZ", "TET3cKO_LZ", "WT_2C", "TET3cKO_2C")
)

df.summary <- DhmR_5hmC_abundance %>%
  group_by(category, Group) %>%
  summarise(
    # sd = sd(abundance, na.rm = TRUE),
    se = sd(abundance) / sqrt(length(abundance)),
    abundance = mean(abundance)
  )
df.summary

library(broom)
library(tidyr)
library(purrr)

df.test <- DhmR_5hmC_abundance %>%
  separate(
    col = "Group",
    sep = "_",
    into = c("Group", "Stage")) %>%
  nest(data = c(-category, -Stage)) %>%
  mutate(data = map(data, ~ wilcox.test(abundance ~ Group, data = .x)),
         data = map(data, tidy)) %>%
  unnest(data) %>%
  mutate(
    qvalue = p.adjust(p.value, "fdr")
  )
df.test

write.table(
  df.test, file = "data/processed/WT_TETcKO.LZY_to_2C.DhmR.sc_5k.abundance_test.tsv",
  col.names = T, row.names = F, sep = "\t", quote = F
)

df.summary %>%
  filter(category!="stable") %>%
  ggplot(aes(category, abundance, fill=Group)) +
  # geom_col(fill = "lightgray", color = "black") +
  geom_errorbar(
    aes(ymin = abundance-se, ymax = abundance+se),
    position = position_dodge(0.9), width=.5)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(y="5hmCpG levels", x="Dynamic DhMRs") +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 0.75)

ggsave(
  "viz/WT_TET3cKO.LZY_to_2C.DhmR.sc_5k.abundance.pdf",
  width = 8, height = 6
)

