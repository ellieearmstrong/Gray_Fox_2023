#Load Libraries
library(tidyverse)
library(data.table)

#set directory
setwd("/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/pi")


#Populations to iterate through
pops = c("SanMiguel","SanClemente","SanNicolas","SantaCatalina","SantaCruz","SantaRosa","Cinereoargenteus")

#Fxn to compute pi and total sites in region
computePi = function(dataFrame, outputColnamePI, outputColnameCount){
  Pi = dataFrame %>%
    summarise(z = sum(PI))
  Count = dim(dataFrame)[1]
  combo = cbind(Pi,Count)
  colnames(combo) = c(outputColnamePI, outputColnameCount)
  return(combo)
}

for (i in pops){
  
  
  #Empty data frame to fill with summary data
  summaryInfo = data.frame()
  
  for (chrom in 1:32){
    
    #Read file with pi per site in 
    sitesPi = fread(file = paste("chrom", chrom, "_", i, ".rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.varInvar.sorted.pi.sites.pi", sep=""))
    
    PIgenomeWide = computePi(sitesPi, "PI_genomeWide_allSites", "CountSites_allSites") 
    
    chromSummary = cbind(chrom,PIgenomeWide)
    
    #Dataframe with all chromosomes
    summaryInfo = bind_rows(summaryInfo, chromSummary)
    
  }

  write.table(summaryInfo, file = paste(i, "_SummaryFile.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
}
