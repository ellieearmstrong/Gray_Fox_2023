library(dplyr)
library(ggplot2)
library(tidyverse)

setwd('~/Documents/Gray_Fox_Project/Gray_Fox_2023/heterozygosity/')

het_foxes <- read.csv('fox_merge_variants_snps_autosome.GM.AN.qualdp.het', header = TRUE, sep = '\t')
metadata <- read.csv('../metadata/fox_metadata.csv', header = TRUE)

fox_palette <- c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", 
                 "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", 
                 "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468")
het_foxes <- het_foxes %>%
  dplyr::rename(Individual = INDV) %>%
  mutate(heterozygosity=((N_SITES-O.HOM.)/2233352514))

het_fox_metadata <- merge(het_foxes, metadata, by="Individual")

#reorder points
het_fox_metadata$Locality <- factor(het_fox_metadata$Locality, levels = c('Golden Gate National Rec Area', 
                                                                          'Santa Monica Mountains','Santa Catalina',
                                                                          'Santa Cruz','San Clemente','San Nicolas',
                                                                          'Santa Rosa','San Miguel'))

#Plot heterozygosity
ggplot(het_fox_metadata, aes(x=Locality, y=heterozygosity, color = Locality)) + 
  scale_color_manual(name = "Locality", values = fox_palette) +
  geom_point(aes(), size = 2) +
  ylab('Observed Heterozygosity') + 
  xlab('') + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=60,hjust = 1))

