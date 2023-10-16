#Load libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen=999)

#Load files
setwd('~/Documents/Gray_Fox_2023/ComputePi/summary')

#read in file for pi
fnames = list.files(pattern = "\\.pi$")
df =  rbindlist(sapply(fnames, read_delim, delim = ',', simplify = FALSE, col_names = c("file", "mean", "sd", "min", "max")))  %>%
  mutate(file = gsub("Input file: /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/pi/","",file),
         CHROM = str_extract(file, "^[^_]+"), #keep everything before first underscore
         CHROM = gsub("chrom", "chr", CHROM),
         Locality = sapply(strsplit(file, "_"), function(x) x[3]),
         Locality = gsub("([a-z])([A-Z])", "\\1 \\2", Locality), 
         RefGenome = str_extract(file, "(?<=_)[^_]+(?=_)"), #keep everything btwn first and second underscore
         file = NULL,
         mean = parse_number(mean),
         sd = parse_number(sd),
         min = parse_number(min),
         max = parse_number(max),
         Locality_ordered=factor(Locality, 
                                 levels = c('San Miguel',
                                            'Santa Rosa',
                                            'Santa Cruz',
                                            'Santa Catalina',
                                            'San Nicolas',
                                            'San Clemente',
                                            'Mainland')))
piGW = df %>%
  group_by(Locality_ordered, RefGenome) %>% #date collected and place
  summarise_at(c("mean"), list(mean = mean)) %>%
  mutate_if(is.numeric, round, digits=6)

#plot 
ggplot(piGW, aes(x=Locality_ordered , y=mean,  fill=RefGenome)) +
  #geom_boxplot(size = 1.5, position = position_dodge(width = 0.9)) +
  #geom_jitter() +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_brewer(palette = "Set1", name = "Reference \n Genome") +
  labs(x= "Locality",y = expression(pi)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        axis.title=element_text(size=32),
        legend.title=element_text(size=28, hjust=0.5),
        legend.text=element_text(size=26),
        strip.text = element_text(size=20),
        legend.position = "right")


