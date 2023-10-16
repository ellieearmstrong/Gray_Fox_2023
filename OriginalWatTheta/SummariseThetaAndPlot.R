#Load libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen=999)

#Load files
setwd('~/Documents/Gray_Fox_2023/OriginalWatTheta/SegSitesPerChrom/')

#Read in meta data and randomly select two samples from each group
set.seed(200)
metadata = read.csv('~/Documents/Gray_Fox_2023/metadata/fox_metadata.csv', header = TRUE) %>%
  mutate(Locality_v2 = gsub("Golden Gate National Rec Area", "Mainland", Locality),
         Locality_v2 = gsub("Santa Monica Mountains", "Mainland", Locality_v2)) %>%
  group_by(Locality_v2) %>%
  slice_sample(n = 2, replace = FALSE) %>%
  ungroup()

# #read in file for theta
# fnames = list.files(pattern = "\\.segSites$")
# for (f in seq_along(fnames)) {
#   
#   df = read_delim(fnames[f], delim = "\t") %>%
#     mutate(INDV = gsub("\\..*", "", INDV)) %>%
#     filter(INDV %in% metadata$Individual) %>%
#     mutate(Locality = metadata$Locality_v2[match(INDV, metadata$Individual)],
#            INDV = NULL) %>%
#     group_by(`#CHROM`, POS, Locality) %>%
#     distinct() %>%
#     ungroup() %>%
#     group_by(`#CHROM`, Locality) %>%
#     count(name = "K_perChrom") %>%
#     mutate(RefGenom = str_extract(fnames[f], "(?<=_)[^_]+(?=_)"))
#   
#   write.table(df, file = "~/Documents/Gray_Fox_2023/OriginalWatTheta/ThetaPerChrom_allRefs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = T)
#   
#   gc()
#   
#   print(f)
# }

#Compute Harmonic mean
compute_constant = function(numberChromsSampled){
  x = seq(1, numberChromsSampled-1, 1)
  constant = sum((1/x))
  return(constant)
}
Newconstant = compute_constant(4)

#get the callable sites
callableSitesPerChrom = read_delim("~/Documents/Gray_Fox_2023/ComputePi/summary/numSitesPerChrom.txt", delim = " ", col_names = c("callableSites", "file")) %>%
  mutate(file = gsub("/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/pi/","",file),
         CHROM = str_extract(file, "^[^_]+"), #keep everything before first underscore
         CHROM = gsub("chrom", "chr", CHROM),
         RefGenome = str_extract(file, "(?<=_)[^_]+(?=_)"), #keep everything btwn first and second underscore
         file = NULL,
         denominator = callableSites*Newconstant)

#Compute genome wide theta
callableSitesGW = data.frame(
  callableSites = c(1480881980,1506020035,1485649587,1483510359),
  RefGenome = c("arcticfox","grayfox","Canfam3.1","Canfam4")
) %>%
  mutate(denominator = callableSites*Newconstant)

#Read in the output file for theta and compute it
summaryPerChrom = read_delim("~/Documents/Gray_Fox_2023/OriginalWatTheta/ThetaPerChrom_allRefs.txt", delim = "\t", col_names = c("CHROM","Locality", "Count", "RefGenome")) %>%
  left_join(callableSitesPerChrom) %>%
  mutate(WattersonsTheta = Count/denominator,
         Locality_ordered=factor(Locality, 
                                 levels = c('San Miguel',
                                            'Santa Rosa',
                                            'Santa Cruz',
                                            'Santa Catalina',
                                            'San Nicolas',
                                            'San Clemente',
                                            'Mainland'))) #mainland vs island
meanPerChrom = summaryPerChrom %>%  
  group_by(Locality_ordered, RefGenome) %>%
  summarise_at(c("WattersonsTheta"), list(mean = mean, sd = sd)) %>%
  mutate_if(is.numeric, round, digits=6)
  

#aggregate genomes wide
meanGW = summaryPerChrom %>%
  group_by(Locality_ordered, RefGenome) %>%
  summarise_at(.vars = c("Count"), sum) %>%
  ungroup() %>%
  left_join(callableSitesGW) %>%
  mutate(WattersonsTheta = Count/denominator) %>%
  mutate_if(is.numeric, round, digits=6)


#Plot results
ggplot(summaryPerChrom, aes(x=Locality_ordered, y=as.numeric(mean), fill=RefGenome)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_brewer(palette = "Set1", name="Reference Genome") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  labs(x= "Locality", y = expression(theta[W])) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 42),
        legend.title = element_text(size = 42),
        legend.text = element_text(size = 40),
        legend.position = "right")

ggtexttable(meanPerChrom, rows = NULL, theme = ttheme("mBlack"))