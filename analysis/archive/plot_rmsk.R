# BiocManager::install("rtracklayer")
library(dplyr)

mouse_TE_age <- read.table("F:/Project/sc5hmC/data/ref/mm10_TE_JCage.tsv", stringsAsFactors = F)

colnames(mouse_TE_age) <- c(
  "chromosome_name", "start_position", "end_position",
  "name", "repFamily", "repClass", "strand",
  "substitution_proportion", "jc_distance"
)

# substitution rate of mouse genome from Waterston et al., 2002
mouse_TE_age$jc_distance <- as.numeric(mouse_TE_age$jc_distance)
mouse_TE_age$mya <- (mouse_TE_age$jc_distance*100)/(4.5*2*100)*1000

# filter repeat annotation file
mouse_TE_age$repFamily <- gsub("\\?","", mouse_TE_age$repFamily)
mouse_TE_age$repClass <- gsub("\\?","", mouse_TE_age$repClass)

mouse_TE_age_filtered %>%
  group_by(repClass, repFamily) %>%
  summarise(mya = mean(mya)) %>%
  arrange(mya)
