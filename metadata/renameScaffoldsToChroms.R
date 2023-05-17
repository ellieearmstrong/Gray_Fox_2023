library(tidyverse)


#read files in
map = read_delim("Foxes/metaData/fox_chrom.map", col_names = F, delim = "\t")
df = read_delim("Foxes/metaData/autosomes_mappability1.0.bed", col_names = F, delim = "\t") %>%
  filter(X1 %in% map$X1) %>%
  mutate(chrom = map$X2[match(X1, map$X1)],
         X1 = NULL) %>%
  na.omit(chrom) %>%
  select(chrom, everything())

gc()

write.table(df, "Foxes/metaData/autosomes_mappability1.0_renameToChroms.bed", col.names = F, row.names = F, quote = F, sep = "\t")
