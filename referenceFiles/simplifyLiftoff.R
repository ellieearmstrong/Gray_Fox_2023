#Load Libraries
library(tidyverse)


#set directory
setwd("~/Documents/Gray_Fox_2023/referenceFiles/")


#Read file with liftoff and chromosomes to rename
chrom_map = read_delim(file = "~/Documents/Gray_Fox_2023/metadata/fox_chrom.map", delim="\t", col_names=c("scaff", "chrom"))

liftoff = read_delim(file = "~/Documents/Gray_Fox_2023/referenceFiles/liftoff_canfam3_grayfox_polished.gff", delim="\t", skip=3, col_names=c("scaff", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>%
  mutate(chrom = chrom_map$chrom[match(scaff, chrom_map$scaff)]) %>%
  na.omit(chrom) %>%
  group_by(chrom, start, end) %>%
  summarise(features=paste(feature, collapse=',')) 
  

write.table(liftoff, "liftoff_canfam3_grayfox_polished_chromsOnlySimplified.out", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

