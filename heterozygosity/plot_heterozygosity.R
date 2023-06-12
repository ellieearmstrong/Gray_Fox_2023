#Load libraries
library(tidyverse)
library(ggpubr)

#Load files
setwd('~/Documents/Gray_Fox_Project/Gray_Fox_2023/heterozygosity/')

het_foxes <- read.delim('fox_merge_variants_snps_autosome.GM.AN.qualdp.het', header = TRUE, sep = '\t')
metadata <- read.csv('../metadata/fox_metadata.csv', header = TRUE)

#Reformat files
het_foxes <- het_foxes %>%
  rename(Individual = INDV) %>%
  mutate(heterozygosity=((N_SITES-O.HOM.)/2233352514))

het_fox_metadata <- merge(het_foxes, metadata, by="Individual") %>%
  mutate(CollectionGeneral = case_when(
    Collection_Date == 1988 ~ "1988", 
    Collection_Date >= 1990 & Collection_Date < 2000 ~ "1990s",
    Collection_Date >= 2000 ~ "2000s"),
    Local_date = paste(Locality, CollectionGeneral),
    Locality_ordered=factor(Locality, 
                               levels = c('San Miguel',
                                          'Santa Rosa',
                                          'Santa Cruz',
                                          'Santa Catalina',
                                          'San Nicolas',
                                          'San Clemente',
                                          'Santa Monica Mountains',
                                          'Golden Gate National Rec Area')),
    hetPer10Kb = heterozygosity*10000)

#Make summary data for bar chart
barPlotSummary = het_fox_metadata %>% 
  group_by(Local_date) %>% 
  summarise_at(.vars = c("hetPer10Kb"), 
               list(mean=mean, sd = sd)) %>%
  mutate(sd = ifelse(is.na(sd), 0, sd),
         Local_date_ordered= factor(Local_date, 
                                    levels = c('Santa Monica Mountains 2000s',
                                               'Golden Gate National Rec Area 1990s',
                                               'San Clemente 1988',
                                               'San Clemente 2000s',
                                               'San Nicolas 1988',
                                               'San Nicolas 2000s',
                                               'Santa Catalina 1988',
                                               'Santa Catalina 2000s',
                                               'Santa Cruz 1988',
                                               'Santa Cruz 2000s',
                                               'Santa Rosa 1988',
                                               'Santa Rosa 2000s',
                                               'San Miguel 1988',
                                               'San Miguel 2000s')),
         CollectionGeneral = het_fox_metadata$CollectionGeneral[match(Local_date, het_fox_metadata$Local_date)])

#Plotting
cbPalette <- c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", 
                 "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", 
                 "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468")

#Plot heterozygosity per base pair
points = ggplot(het_fox_metadata, aes(x=Locality_ordered, y=heterozygosity, color = Locality)) + 
  scale_color_manual(name = "Locality", values = cbPalette) +
  geom_point(size = 8) +
  ylab('Observed Heterozygosity\n(per base pair)') + 
  xlab(NULL) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 30, angle = 90, vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(size = 30), 
        plot.title=element_text(size=32, face = "bold", hjust=0.5), 
        axis.title=element_text(size=32),
        legend.position = "none")

#Plot heterozygosity bar
bar = ggplot(data = barPlotSummary, aes(fill = CollectionGeneral)) + 
  scale_fill_brewer(name = "Date", palette = "Dark2") +
  geom_bar(aes(x=Local_date_ordered, y=mean), stat="identity", size = 2, position = dodge) +
  geom_errorbar(aes(x=Local_date_ordered, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.3) +
  coord_flip() +
  ylab('Autosomal Heterozygosity\n(per 10kb)') + 
  xlab(NULL) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        plot.title=element_text(size=32, face = "bold", hjust=0.5), 
        axis.title=element_text(size=32),
        legend.position = "none")

#put plots together
ggarrange(points, bar, nrow = 1)

