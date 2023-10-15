#Load libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen=999)

#Load files
setwd('~/Documents/Gray_Fox_2023/heterozygosity/')
metadata = read.csv('~/Documents/Gray_Fox_2023/metadata/fox_metadata_extended.csv', header = TRUE) %>% 
  select(-c(heterozygosity, F, hetPer10Kb))
fnames = list.files(pattern = "\\.vcf.gz.het$")

#Generate data frame
##create columns with fileName, population, and compute pi 
df = rbindlist(sapply(fnames, read_delim, delim = '\t', simplify = FALSE), use.names = TRUE, idcol = "Ref") %>%
  mutate(Ref = gsub("_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz.het","", Ref),
         INDV = gsub("\\..*", "", INDV),
         heterozygosity = case_when(
           Ref == "Canfam4" ~ (N_SITES-`O(HOM)`)/1483510359,
           Ref == "Canfam3.1" ~ (N_SITES-`O(HOM)`)/1485649587,
           Ref == "arcticfox" ~ (N_SITES-`O(HOM)`)/1480881980,
           Ref == "grayfox" ~ (N_SITES-`O(HOM)`)/1506020035,
         )) %>%
  left_join(metadata) %>%
  mutate(hetPer10Kb = heterozygosity*10000,
         Locality_ordered=factor(Locality, 
                               levels = c('San Miguel',
                                          'Santa Rosa',
                                          'Santa Cruz',
                                          'Santa Catalina',
                                          'San Nicolas',
                                          'San Clemente',
                                          'Santa Monica Mountains',
                                          'Golden Gate National Rec Area')))

#Make summary data for bar chart
barPlotSummary = df %>% 
  group_by(Ref, Local_date) %>% 
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
         CollectionGeneral = df$CollectionGeneral[match(Local_date, df$Local_date)])

#compute mean a couple different ways
summaryDF = df %>%
  group_by(Ref, Type) %>% #mainland vs island
  summarise_at(c("heterozygosity", "hetPer10Kb"), list(mean = mean, sd = sd))

summaryDF_date = df %>%
  group_by(Ref, Type, CollectionGeneral) %>% #date collected and place
  summarise_at(c("heterozygosity", "hetPer10Kb"), list(mean = mean, sd = sd))

#Plotting
cbPalette = c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468")

#Plot heterozygosity per base pair across all references
#pdf(file = "Heterozygosity_Indv_barplot.pdf", width = 20, height = 10)
ggplot(df, aes(x=INDV, y=heterozygosity, fill = Ref)) + 
  #facet_wrap(INDV ~ .) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_brewer(name = "Reference Genome", palette = "Set1") +
  ylab('Observed Heterozygosity\n(per base pair)') + 
  xlab(NULL) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 30), 
        axis.title=element_text(size=32),
        legend.title=element_text(size=28, hjust=0.5),
        legend.text=element_text(size=26),
        strip.text = element_text(size=20),
        legend.position = "right")
#dev.off()

#Bar plot of heterozygsity of island foxes based on collection time
#pdf(file = "Heterozygosity_Island_barplot.pdf", width = 20, height = 10)
ggplot(summaryDF_date %>% filter(Type == "Island"), aes(x=CollectionGeneral, y=as.numeric(heterozygosity_mean), fill=Ref)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_brewer(palette = "Set1", name="Reference\nGenome") +
  geom_errorbar(aes(ymin=heterozygosity_mean-heterozygosity_sd, ymax=heterozygosity_mean+heterozygosity_sd), width=.2,
                position=position_dodge(.9)) +
  labs(x = "Collection Time", y="Observed Heterozygosity\n(per base pair") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        legend.title = element_text(size = 42),
        legend.text = element_text(size = 40),
        legend.position = "right")
#dev.off()

#Bar plot heterozygosity in mainland versus island
#pdf(file = "Heterozygosity_IslandMainland_barplot.pdf", width = 20, height = 10)

ggplot(summaryDF, aes(x=Type, y=as.numeric(heterozygosity_mean), fill=Ref)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_brewer(palette = "Set1", name="Reference Genome") +
  geom_errorbar(aes(ymin=heterozygosity_mean-heterozygosity_sd, ymax=heterozygosity_mean+heterozygosity_sd), width=.2,
                position=position_dodge(.9)) +
  labs(x = "Location", y="Observed Heterozygosity\n(per base pair") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        legend.title = element_text(size = 42),
        legend.text = element_text(size = 40),
        legend.position = c(0.3, 0.8))

#dev.off()

#Plot heterozygosity per base pair
points = ggplot(df, aes(x=Locality_ordered, y=heterozygosity, color = Locality)) + 
  facet_wrap(Ref ~ .) +
  scale_color_manual(name = "Locality", values = cbPalette) +
  geom_point(size = 8) +
  ylab('Observed Heterozygosity\n(per base pair)') + 
  xlab(NULL) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 30), 
        axis.title=element_text(size=32),
        legend.title=element_text(size=28, hjust=0.5),
        legend.text=element_text(size=26),
        strip.text = element_text(size=20),
        legend.position = "right")
print(points)
#pdf(file = "Heterozygosity_perLocality_perRef.pdf", width = 20, height = 10)
#print(points)
#dev.off()

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

