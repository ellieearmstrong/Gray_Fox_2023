#Load libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen=999)

#Load files
setwd('~/Documents/Gray_Fox_2023/OriginalWatTheta/SegSitesPerChrom/')
fnames = list.files(pattern = "\\.segSites$")

#Read in meta data and randomly select two samples from each group
set.seed(200)
metadata = read.csv('~/Documents/Gray_Fox_2023/metadata/fox_metadata.csv', header = TRUE) %>%
  mutate(Locality_v2 = gsub("Golden Gate National Rec Area", "Mainland", Locality),
         Locality_v2 = gsub("Santa Monica Mountains", "Mainland", Locality_v2)) %>%
  group_by(Locality_v2) %>%
  slice_sample(n = 2, replace = FALSE) %>%
  ungroup()

# #read in file for theta
# for (f in seq_along(fnames)) {
#   
#   df = read_delim(fnames[f], delim = "\t") %>%
#     mutate(INDV = gsub("\\..*", "", INDV)) %>%
#     filter(INDV %in% metadata$Individual) %>%
#     mutate(Locality = metadata$Locality_v2[match(INDV, metadata$Individual)],
#            INDV = NULL) %>%
#     group_by(`#CHROM`, POS, Locality) %>%
#     distinct() %>%
#     ungroup() %>%
#     group_by(`#CHROM`, Locality) %>%
#     count(name = "K_perChrom") %>%
#     mutate(RefGenom = str_extract(fnames[f], "(?<=_)[^_]+(?=_)"))
#   
#   write.table(df, file = "~/Documents/Gray_Fox_2023/OriginalWatTheta/ThetaPerChrom_allRefs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = T)
#   
#   gc()
#   
#   print(f)
# }

#Read in the output file
summary = read_delim("~/Documents/Gray_Fox_2023/OriginalWatTheta/ThetaPerChrom_allRefs.txt", delim = "\t", col_names = c("CHROM","Locality", "Count", "RefGenome"))

