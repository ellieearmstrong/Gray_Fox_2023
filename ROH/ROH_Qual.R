#Load libraries
library(data.table)
library(tidyverse)
library(mgsub)

#Load files
setwd("~/Documents/Gray_Fox_2023/ROH/overlapCallable/")
popDF = read.csv("~/Documents/Gray_Fox_2023/metadata/fox_metadata.csv")
chroms = read.delim("~/Documents/Gray_Fox_2023/metadata/fox_chrom_withLengths.map")

#Load vcftools outputs
fnames = list.files(pattern = "\\ROH.bed$")
df = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE) #Read all the files into a data frame 
colnames(df) = c("CHROM", "start", "end", "type" ,"indiv", "numCallableSitesOverlap", "numROHBasesOverlaped", "rohLength", "propCoverage")
#Remove header lines and convert character columns to numeric

filterROH = df %>%
  select("CHROM", "start","type" , "end", "indiv", "rohLength", "propCoverage") %>%
  filter(propCoverage >= .60 & propCoverage <= .95 & rohLength >= 10000) %>%
  mutate(Population = popDF$Locality[match(indiv, popDF$Individual)])
  

#Alternative way to QC ROH Segments
##Will remove ROH if the proportion of ROH segment covered by snps is not within a SD of the mean
#z = data.table(df)
#z[,ToKeep := abs(df$propCoverage - mean(df$propCoverage)) < sd(df$propCoverage)][ToKeep  == TRUE] #create variable that identifies what to drop or keep 
#filterROH = subset(z, z$ToKeep == "TRUE")


#write.table(filterROH, file="~/Documents/Gray_Fox_2023/ROH/FinalROHQCd_min10KbminPropCov60per_allIndivs_May2023_garlic.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Plot the distribution of coverage of ROHs post-filtering
postFilter = ggplot() + 
  geom_density(data = filterROH, aes(x=propCoverage)) + 
  labs(x="Proportion of ROH Covered") + 
  ggtitle("Post-Filtering") + 
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24))
#print(postFilter)

#remove anything less than 1Mb and approximate FROH
genomeLength = sum(chroms$length)

IndividualLevel = filterROH %>%
  filter(type == "C") %>%
  group_by(indiv) %>% 
  summarise(TotalROH = sum(rohLength)) %>%
  mutate(TotalROHMb = TotalROH/10^6, 
         FROH = TotalROH/genomeLength,
         Population = popDF$Locality[match(indiv, popDF$Individual)])



#Plot FROH with ROH greater than 1Mb
#prep data for plotting
cbPalette <- c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", 
               "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", 
               "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468", "Cinereoargenteus" = "gray20")

FROHplot = ggplot(IndividualLevel, aes(x=Population, y=FROH, colour=Population)) + 
  scale_colour_manual(name = "Population", values = cbPalette) + 
  geom_boxplot(size=1) + 
  geom_point(show.legend=F) + 
  theme_bw() + 
  labs(x = "Population", y=expression(F[ROH])) + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

