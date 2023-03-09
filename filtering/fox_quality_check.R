library(dplyr)
library(ggplot2)
#setwd
setwd('~/Documents/Gray_Fox_Project/Gray_Fox_2023/')

quality_foxes <- read.table('~/Documents/Gray_Fox_Project/Gray_Fox_2023/grayfox_quality.tab')

quality_foxes <- quality_foxes %>%
  rename('chrom' = 'V1', 'position' = 'V2', 'ref' = 'V3', 'alt' = 'V4',
         'AN' = 'V5', 'DP' = 'V6')


mean(quality_foxes$AN)
mean(quality_foxes$DP)
summary(quality_foxes$DP)

a <- ggplot(quality_foxes, aes(DP)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

a + theme_light() + xlim(0, 1000)
