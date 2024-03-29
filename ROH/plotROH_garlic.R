#Load libraries and set working directory
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("~/Documents/Tigers/ROH/")

#Load roh and FSNP and meta data
popsDF = read_delim("~/Documents/Tigers/IndivFiles/TableX-SampleDetails.txt")

roh = read_delim("~/Documents/Tigers/ROH/TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt", delim = "\t") %>%
  left_join(popsDF, by = c("INDV" = "Sample")) %>%
  mutate(Subspecies2 = factor(Subspecies_GroupID_Corrected, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')),
         TYPE2 = paste("Type",TYPE))

rohLengthsClass = roh %>%
  group_by(INDV, TYPE) %>%
  summarise_at(c("AUTO_LEN"), sum) %>%
  ungroup() %>%
  left_join(popsDF,by = c("INDV" = "Individual")) %>%
  mutate(Subspecies2 = factor(Subspecies2, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')),
         TYPE2 = paste("Type",TYPE))

FROH = rohLengthsClass %>%
  filter(TYPE == "C") %>%
  mutate(Froh = AUTO_LEN/2174711735) 

summaryTable = roh %>%
  mutate(TYPE2 = paste("Type",TYPE)) %>%
  arrange(TYPE2, Subspecies2) %>%
  group_by(TYPE2, Subspecies2) %>%
  summarise(mean = mean(AUTO_LEN)/10^3, 
            min = min(AUTO_LEN)/10^3, 
            max = max(AUTO_LEN)/10^3) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  ungroup() %>%
  rename("ROH class" = TYPE2, "Subspecies" = "Subspecies2") 

ggtexttable(summaryTable, rows = NULL, theme = ttheme("mBlackWhite"))

##Plot data
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Generic"="gray25")#palette
cbPalette_expanded = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Generic"="gray25", "Generic-Orange" = "#CC79A7", "Generic-SnowWhite" = "#867BCF", "Generic-Golden"="darkseagreen3", "Generic-White"="cornflowerblue")#palette


plotROHs = ggplot(rohLengthsClass, aes(x=Subspecies2, y=AUTO_LEN/10^6, fill=Subspecies2)) + 
  geom_violin() +
  geom_jitter(height = 0, width = 0.1) +
  facet_wrap(~TYPE2) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) +
  labs(x="Subspecies", y="Length of genome in ROH (Mb)") +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))  
print(plotROHs)

plotFROH = ggplot(FROH, aes(x=Subspecies2, y=Froh)) +
  geom_violin(aes(fill=Subspecies2)) + 
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) + 
  labs(x="Subspecies", y=expression(F[ROH])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14 )) 

print(plotFROH)

###plot inbreeding measures
y = all_nodups_full_highCov %>%
  select(Individual, FSNP, Subspecies2)
x = FROH %>%
  select(INDV, Froh, Subspecies2) 
colnames(x) = c("INDV", "Value","Subspecies2" ,"TYPE")
colnames(y) = c("INDV", "Value","Subspecies2" ,"TYPE")
z = rbind.data.frame(x, y)
z$facets = factor(z$TYPE, 
                  labels = c("F[ROH]", "F[SNP]"))
inbreeding = ggplot(z, aes(x=Subspecies2, y=Value, fill=Subspecies2)) + 
  geom_violin() +
  geom_jitter(height = 0, width = 0.1) +
  facet_wrap(~facets, labeller = label_parsed) +
  coord_flip() +
  scale_y_continuous(breaks=seq(-0.5,1,0.25)) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) +
  labs(x="Subspecies", y="Consanguinity coefficient value") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))  

plot = ggarrange(plotROHs + xlab(NULL), 
                 inbreeding + xlab(NULL), 
                 nrow = 2, 
                 common.legend = TRUE, 
                 legend = "none")

annotate_figure(plot, 
                left = text_grob("Subspecies", color = "black", face = "bold", size = 18, rot = 90))
