#Load libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen=999)

####Function to make a folded SFS###
#THIS IS BERNARD's ORIGINAL >> need to add 0 to front and end to stand in for monomorphic (don't need to do this if it came from dadi)
fold <- function(SFSCountCol, n, norm=TRUE){
  if (length(SFSCountCol) < (n+1)){
    data = c(SFSCountCol, rep(0,(n+1)-length(SFSCountCol)))
  }
  data = SFSCountCol[2:n] # adds together sfs and backwards sfs
  data_fold = data + rev(data) # takes first half of that added together sfs (but not middle entry)
  data_fold = data_fold[1:(n/2-1)]# adds middle entry that didn't have anything added to the end
  data_fold = c(data_fold,data[(n/2)])# truncates and sums up remaining fields if desired (not needed here)
  #data_trunc = c(data_fold[1:(trunc-1)],sum(data_fold[trunc:length(data_fold)]))
  #if (norm){
  #  data_trunc = data_trunc/sum(data_trunc)
  #}
  #return(data_trunc)
  return(tibble(data_fold))
}

####Function to compute harmonic mean of chroms for watterson's theta #####
compute_constant = function(numberChromsSampled){
  x = seq(1, numberChromsSampled-1, 1)
  constant = sum((1/x))
  return(constant)
}

####Load files
setwd('~/Documents/Gray_Fox_2023/SFS')
fnames = list.files(pattern ="*_filtered.renameChroms.Mainland.drop295.*")

#Generate data frame
df = rbindlist(sapply(fnames, read.delim, sep="\t", simplify = FALSE), use.names = TRUE, idcol = "filename") %>%
  mutate(ref = gsub("_.*$","",filename),
         pop = gsub(".*_(East|West|Hybrid)_.*", "\\1", filename),
         filename = NULL)

####SFS data 
PlottingFolded_N12 = df %>%
  group_by(ref, pop) %>%
  group_modify(~ fold(.x$countRef, n = 24)) %>% #this n = num chroms
  mutate(proportion = data_fold / sum(data_fold)) %>%
  ungroup() %>%
  mutate(bin = rep(seq(1:12), length.out = n()),
         ref = factor(ref, levels = c('arcticfox', 'grayfox', 'Canfam4', 'Canfam3.1', 'graywolf')))

####Plot SFS
ggplot(PlottingFolded_N12, aes(y=proportion, x=bin, fill=ref)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  facet_grid(~pop) +
  scale_x_continuous(breaks=1:12) + 
  ylim(0,0.35)+
  scale_fill_brewer(palette = "Set1", name = "Reference genome") + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) 

#write.table(PlottingFolded_N12, file = "SFS_allComps.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

####Compute Watterson's Theta from SFS and original data
CallableSites = df %>%
  group_by(ref, pop) %>%
  summarise(TotalSites = sum(countRef)) %>%
  ungroup() 

SegSites = PlottingFolded_N12 %>%
  group_by(ref, pop) %>%
  summarise(SegSites = sum(data_fold)) %>%
  ungroup() 

TotalSites = left_join(CallableSites, SegSites)

####Plot SFS
PlotTheta = TotalSites %>%
  group_by(ref, pop) %>%
  mutate(theta = SegSites/(TotalSites*compute_constant(24))) %>%
  ungroup() %>%
  mutate(ref = factor(ref, levels = c('arcticfox', 'grayfox', 'Canfam4', 'Canfam3.1', 'graywolf')))

ggplot(PlotTheta, aes(y=theta, x=ref, fill=ref)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  facet_grid(~pop) +
  ylim(0,0.2)+
  scale_fill_brewer(palette = "Set1", name = "Reference genome") + 
  labs(x = expression(theta[W]), 
       y= "Reference genome") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) 

#write.table(PlotTheta, file = "WatTheta_allComps.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
