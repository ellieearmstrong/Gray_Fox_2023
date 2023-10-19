#Load libraries
library(data.table)
library(tidyverse)
library(mgsub)

#Load files
setwd("~/Documents/Gray_Fox_2023/ROH/overlapCallable/")
popDF = read.csv("~/Gray_Fox_2023/metadata/fox_metadata_extended.csv")
chroms = read.delim("~/Documents/Gray_Fox_2023/metadata/fox_chrom_withLengths.map")

# #Load vcftools outputs
# fnames = list.files(pattern = "\\ROH.bed$")
# df = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE) #Read all the files into a data frame 
# colnames(df) = c("CHROM", "start", "end", "INDV", "numCallableSitesOverlap", "numROHBasesOverlaped", "rohLength", "propCoverage")

# #filter on proportion coverage 
#  filterROHCoverage = df %>%
#    select("CHROM", "start", "end", "INDV", "rohLength", "propCoverage") %>%
#    filter(propCoverage >= .60 & rohLength >= 10000) %>% #60 % coverage and at least 10kb
#    mutate(Ref = gsub('.*\\.', '', INDV),
#           Ref = gsub("1", "Canfam3.1", Ref),
#           INDV = gsub("\\..*", "", INDV)) %>%
#    left_join(allfiles) %>%
#    left_join(popDF) 
#  
# #find shared ROH, mark to drop if more than two individuals have them  
#  sharedROH = filterROHCoverage %>% 
#    group_by(CHROM, start, end, Ref) %>% 
#    count() %>% 
#    mutate(keep = ifelse(n < 3, "keep", "toss")) 
#  
# #Remove shared stuff
#  filterROH = filterROHCoverage %>%
#    left_join(sharedROH) %>%
#    filter(keep == "keep") %>%
#    mutate(keep = NULL,
#           Locality_ordered=factor(Locality, 
#                                   levels = c('San Miguel',
#                                              'Santa Rosa',
#                                              'Santa Cruz',
#                                              'Santa Catalina',
#                                              'San Nicolas',
#                                              'San Clemente',
#                                              'Santa Monica Mountains',
#                                              'Golden Gate National Rec Area')))

#Read in file with filtered ROH 
 filterROH = read_delim("FinalROHQCd_min10KbminPropCov60per_allIndivs_May2023_garlic.bed", delim = "\t") %>%
   mutate(keep = NULL,
          Locality_ordered=factor(Locality, 
                                  levels = c('San Miguel',
                                             'Santa Rosa',
                                             'Santa Cruz',
                                             'Santa Catalina',
                                             'San Nicolas',
                                             'San Clemente',
                                             'Santa Monica Mountains',
                                             'Golden Gate National Rec Area')))


#write.table(filterROH, file="~/Gray_Fox_2023/ROH/FinalROHQCd_PropCovgrThan60per_allIndivs_Oct2023_garlic.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

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


#Compute FROH with the total ungapped length
FROH = filterROH %>%
  filter(ROHCategory == "C") %>%
  group_by(Ref, INDV) %>% 
  summarise(TotalROH = sum(rohLength)) %>%
  mutate(FROH = case_when(
           Ref == "Canfam4" ~ TotalROH/2482000080,
           Ref == "Canfam3.1" ~ TotalROH/2392715236,
           Ref == "arcticfox" ~ TotalROH/2345550353,
           Ref == "grayfox" ~ TotalROH/2658766243)) %>%
  left_join(popDF) %>%
  mutate(Locality_ordered=factor(Locality, 
                                 levels = c('San Miguel',
                                            'Santa Rosa',
                                            'Santa Cruz',
                                            'Santa Catalina',
                                            'San Nicolas',
                                            'San Clemente',
                                            'Santa Monica Mountains',
                                            'Golden Gate National Rec Area')))



#Plotting
#color palette
cbPalette <- c("Golden Gate National Rec Area" = "#FBD05C", "Santa Monica Mountains" = "#B2A6CC", 
               "Santa Catalina" ="#8FFFF5", "Santa Cruz" = "#4C8C42", "San Clemente" = "#3854A6", 
               "San Nicolas"= "#F25C05", "Santa Rosa" = "#B20650", "San Miguel" = "#6F4468", "Cinereoargenteus" = "gray20")

FROHplot = ggplot(FROH, aes(x=Locality_ordered, y=FROH, color = Locality)) + 
  facet_wrap(Ref ~ .) +
  geom_point(size = 8) + 
  theme_bw() + 
  labs(x = "Locality", y=expression(F[ROH])) +
  scale_color_manual(name = "Locality", values = cbPalette) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        strip.text = element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18),
        legend.position = "right")
print(FROHplot)


plotROH = ggplot(filterROH, aes(x=Locality_ordered, y=rohLength/10^6, colour = Locality)) + 
  facet_grid(Ref ~ ROHCategory, scales = "free_x") +
  geom_boxplot() + 
  geom_point(color = "black", show.legend = F) +
  coord_flip() +
  theme_bw() + 
  labs(y = "Locality", x="ROH Length (Mb)") +
  scale_colour_manual(name = "Locality", values = cbPalette) +
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        strip.text = element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18),
        legend.position = "none")
print(plotROH)