library(tidyverse)
library(data.table)

setwd("~/proj/Jazlyn/grayfox/Analyses/pi/output_pi_50Kb")

refs = unique(str_split_i(list.files(pattern = "*filtered.renameChroms.Mainland.drop295.*"), "_", 1))

for (i in refs) {
  
  fnames = list.files(pattern = paste0("^", i, "_filtered.renameChroms.Mainland.drop295.*"))
  
  #Generate data frame
  df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE) %>%
    na.omit() %>%
    filter(count_missing < 7000) %>% #at least 200 snps in the 50kb window
    mutate(chromosome = as.numeric(gsub("chr","",chromosome)))
  
  avg_windowed_pi_perChrom = df %>%
    group_by(pop,chromosome) %>%
    summarise_at(.vars = c("avg_pi"), list(mean = mean, sd = sd))
  
  write.table(avg_windowed_pi_perChrom, file = paste0("pi_", i, "_allComps_allChroms.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  avg_windowed_pi = df %>%
    group_by(pop) %>%
    summarise_at(.vars = c("avg_pi"), list(mean = mean, sd = sd))
  
  write.table(avg_windowed_pi, file = paste0("pi_", i, "_allComps.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  cat(sprintf("done %s", i),"\n")
}


