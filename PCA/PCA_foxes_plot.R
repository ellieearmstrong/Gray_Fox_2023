library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)

setwd('~/Documents/Gray_Fox_Project/Gray_Fox_2023/PCA/')

eigenvec_fox <- read.table('fox_merge_variants_snps_autosome.GM.AN.qualdp.varid.eigenvec',head=F)
eigenval_fox <- read.table('fox_merge_variants_snps_autosome.GM.AN.qualdp.varid.eigenval',head=F)
metadata <- read.csv('../metadata/fox_metadata.csv', header = TRUE)

fox_palette <- c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", 
                 "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", 
                 "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468")


eigenval_fox$V2 <- (eigenval_fox$V1/16)*100

eigenvec_fox <- eigenvec_fox %>% 
  rename(Individual=V2)
eigenvec_fox<- merge(eigenvec_fox, metadata,by="Individual")

PC1 <- ggplot(eigenvec_fox, aes(x=V3, y=V4, color=Locality)) +
  scale_color_manual(name="Locality", values=fox_palette) +
  geom_point(size = 2) +
  xlab("PC1 (30.7%)") + ylab("PC2 (19.5%)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme_bw()

PC2 <- ggplot(eigenvec_fox, aes(x=V3, y=V5, color=Locality)) +
  scale_color_manual(name="Locality", values=fox_palette) +
  geom_point() +
  xlab("PC1 (30.7%)") + ylab("PC3 (17.7%)") +
  theme_bw() +
  theme(legend.position = "none") 

PC3 <- ggplot(eigenvec_fox, aes(x=V4, y=V5, color=Locality)) +
  scale_color_manual(name="Locality", values=fox_palette) +
  geom_point() +
  xlab("PC2 (19.5%)") + ylab("PC3 (17.7%)") +
  theme_bw() +
  theme(legend.position = "none")

PC1 | 
  (PC2/PC3)  

