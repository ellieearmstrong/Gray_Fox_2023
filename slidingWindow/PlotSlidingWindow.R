#Load Libraries
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(ggpubr)

#Load files
setwd("~/Documents/Gray_Fox_2023/slidingWindow/1000000win_1000000step/")

#Function will annotate files so that have coverage percentages 
#Window size is 1 million so we want coverage percentage based on 10^6
annotateDF = function(inputDF){
  inputDF %>%
    filter(sites_passing != 0) %>%
    mutate(chromo = as.numeric(as.character(gsub("chr", "", chromo))), 
           window_end = window_start + 10^6, 
           TwentyPerCov = ifelse(sites_total >= 2e+05, "1", "0"),
           FiftyPerCov = ifelse(sites_total >= 5e+05, "1", "0"), 
           SixtyPerCov = ifelse(sites_total >= 6e+05, "1", "0")) 
}

#This function will make a new data frame and adjust windows
ProcessDF = function(inputDF){ inputDF %>% 
    
    # Compute chromosome size
    group_by(chromo) %>% 
    summarise(chr_len=max(window_start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(inputDF, ., by=c("chromo"="chromo")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chromo, window_start) %>%
    mutate( newWinStart=window_start+tot) }

#This function implements the previous two functions to create our final data frame for plotting

makeFinalFiles = function(filesToJoin){
  Fnames = lapply(Sys.glob(filesToJoin), read.delim) #grab all chroms
  df = rbindlist(Fnames) #merge all chroms
  dfAnnot = annotateDF(df) %>% 
    filter(FiftyPerCov == 1) #use only windows with at least 50 percent coverage if you increase TM doesn't have very many windows
  processedDF = ProcessDF(dfAnnot)
  return(processedDF)
}

#Generate dataframes from all sliding windows
fox_Annotated = makeFinalFiles("chrom*_het_1000000win_1000000step.txt")

#Function to plot data for each individual
plotFunction = function(dataFrame, indiv, color1, color2) {
  #Generate x axis with any one of the data frames 
  axisdf = dataFrame %>% 
    group_by(chromo) %>% 
    summarize(center=(max(newWinStart) + min(newWinStart) ) / 2 )
  #Now plot with the axis
  indivHet = ggplot() + 
    geom_bar(data = dataFrame, 
             aes(x=newWinStart, y=(indiv/sites_total)*1000, 
                 color=as.factor(chromo)),
             stat = "identity", 
             lwd=0.5) +
    scale_color_manual(values = rep(c(color1, color2), 32)) +
    #custom X axis:
    scale_x_continuous(label = axisdf$chromo, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0)) + # remove space between plot area and x axis
    labs(x = "Chromosome", y = "Heterozygosity (per kb)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          panel.border = element_blank(), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.position = "none")
  return(indivHet)
}

#Plot Santa Catalina
hets_SRR7458271 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458271, "gray50", "#8FFFF5")
hets_SCA16F = plotFunction(fox_Annotated, fox_Annotated$hets_SCA16F, "gray50", "#8FFFF5")
#Plot San Clemente
hets_SCLV4F = plotFunction(fox_Annotated, fox_Annotated$hets_SCLV4F, "gray50", "#3854A6")
hets_SRR7458268 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458268, "gray50", "#3854A6")
#Plot Santa Cruz
hets_SCZ05M = plotFunction(fox_Annotated, fox_Annotated$hets_SCZ05M, "gray50", "#4C8C42")
hets_SRR7458269 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458269, "gray50", "#4C8C42")
#Plot San Miguel
hets_SMI15F = plotFunction(fox_Annotated, fox_Annotated$hets_SMI15F, "gray50", "#6F4468")
hets_SRR7458266 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458266, "gray50", "#6F4468")
#Plot Santa Rosa
hets_SRO40F = plotFunction(fox_Annotated, fox_Annotated$hets_SRO40F, "gray50", "#B20650")
hets_SRR7458265 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458265, "gray50", "#B20650")
hets_SRR7458264 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458264, "gray50", "#B20650")
#Plot Santa Monica Mnts
hets_GFO41F = plotFunction(fox_Annotated, fox_Annotated$hets_GFO41F, "gray50", "#B2A6CC")
#Plot San Nicolas
hets_SNI05F = plotFunction(fox_Annotated, fox_Annotated$hets_SNI05F, "gray50", "#F25C05")
hets_SNI41F = plotFunction(fox_Annotated, fox_Annotated$hets_SNI41F, "gray50", "#F25C05")
hets_SRR7458267 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458267, "gray50", "#F25C05")
#Plot Golden Gate
hets_SRR7458270 = plotFunction(fox_Annotated, fox_Annotated$hets_SRR7458270, "gray50", "goldenrod")

ggarrange(hets_SRR7458266,hets_SMI15F,hets_SRR7458264,hets_SRR7458265,hets_SRO40F,hets_SRR7458267,hets_SNI41F, hets_SNI05F, hets_SRR7458268,hets_SCLV4F,hets_SRR7458269,hets_SCZ05M, hets_SCA16F, hets_SRR7458271,hets_GFO41F, hets_SRR7458270, nrow = 4, ncol = 4, align='hv')
