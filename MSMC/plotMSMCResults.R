#Set Libraries and working directory
library(tidyverse)
library(scales)
library(data.table)
library(mgsub)
setwd("~/Gray_Fox_2023/MSMC")
fnames = list.files(pattern = "\\.final.txt$") #load up files
#load meta data
popDF = read.csv("~/Gray_Fox_2023/metadata/fox_metadata.csv")

####MSMC outputs times and rates scaled by the mutation rate per basepair per generation. First, scaled times are given in units of the per-generation mutation rate. This means that in order to convert scaled times to generations, divide them by the mutation rate. In humans, we used mu=1.25e-8 per basepair per generation.To convert generations into years, multiply by the generation time, for which we used 30 years.
####To get population sizes out of coalescence rates, first take the inverse of the coalescence rate, scaledPopSize = 1 / lambda. Then divide this scaled population size by 2*mu (yes, this factor 2 is different from the time scaling, sorry).

mu=2.0e-8#Tanya Autosomal Mutation Rate for Neutral cm >0.4 and 3years/gen
gen=1 #(years/gen)

results25 = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName") %>%
  mutate(Ne = (1/lambda)/(2*mu), #note the factor of 2! (not in time scaling) confirmed correct: https://github.com/stschiff/msmc-tools/blob/master/plot_utils.py
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/3,
         RightYears = gen*(right_time_boundary/mu),
         Individual = mgsub(FileName, pattern = c(".*_", ".final.txt"),replacement = c("","")),
         Population = popDF$Locality[match(Individual, popDF$Individual)]) 

#Get parameters for priors 
cbPalette <- c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", 
               "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", 
               "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468", "Cinereoargenteus" = "gray20")

ggplot(results25, aes(x=LeftYears,y=Ne, colour=Population))+
  geom_step(stat="identity")+
  scale_colour_manual(name = "Population", values = cbPalette) + 
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("PSMC Results\nmu= ",mu,"\ngeneration time = ",gen," yrs/gen",sep=""))+
  xlab("Years Ago")+
  ylab("Ne")+
  scale_y_log10(labels=comma,breaks=c(1000,10000,100,000,1000000, 10000000, 100000000, 1000000000, 10000000000))+
  scale_x_log10(breaks=c(100,1000,10000,100,000,1000000, 10000000, 100000000, 1000000000, 10000000000),labels=comma)+
  theme(legend.position="bottom")


ggplot(results25 %>% filter(Population == "Golden Gate National Rec Area",
                            LeftYears < 10^6), aes(x=LeftYears,y=Ne, colour=Population))+
  geom_step(stat="identity")+
  scale_colour_manual(name = "Population", values = cbPalette) + 
  scale_y_log10(labels=comma,breaks=c(1000,10000,100,000,1000000, 10000000, 100000000, 1000000000, 10000000000))+
  scale_x_log10(breaks=c(100,1000,10000,100,000,1000000, 10000000, 100000000, 1000000000, 10000000000),labels=comma) + 
  labs(title = paste("PSMC Results\nmu= ",mu,"\ngeneration time = ",gen," yrs/gen",sep=""),
       x = "Years ago", 
       y=expression(N[E])) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=20, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_blank(), 
        legend.text=element_text(size=18))
