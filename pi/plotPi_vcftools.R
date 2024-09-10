#Load libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen=999)

#read file in
setwd("~/Documents/Gray_Fox_2023/pi") #set working directory

df = read_delim("piPerChrom_allComps.pi", delim = "\t", col_names = c("pi", "info")) %>%
  separate_wider_delim(info,"_", names = c("ref", "chrom", "pop", "N")) %>%
  mutate(N=NULL)

#compute average
avg_pi = df %>%
  group_by(ref, pop) %>%
  summarise(mean = mean(pi, na.rm = TRUE), 
            sd = sd(pi, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ref = factor(ref, levels = c('arcticfox', 'grayfox', 'Canfam4', 'Canfam3.1', 'graywolf')))


#plot
ggplot(avg_pi %>% filter(pop != "Hybrid"), aes(y=mean, x=ref, fill=ref)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_grid(~pop) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(from=0,to=0.16,by=0.04)) +
  scale_fill_brewer(palette = "Set1", name = "Reference genome") + 
  labs(x = "Reference genome", 
       y= expression(pi)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) 
