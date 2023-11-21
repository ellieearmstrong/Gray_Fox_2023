library(ggplot2)
library(patchwork)
library(dplyr)
greyfox <- read.csv('~/Documents/Gray_Fox_Project/GrayFox_RefBiasProject/MappingStats_Plot.csv')


depth <- greyfox %>% 
  filter(Type == 'Depth') 

auto_depth <- greyfox %>% 
  filter(Type == 'AutosomeDepth') %>%
  filter(Dataset == 'Robinson')

auto_depth_gray <- greyfox %>% 
  filter(Type == 'AutosomeDepth') %>%
  filter(Dataset == 'Sacks')

mapped <- greyfox %>% 
  filter(Type == 'Mapping') %>%
  filter(Dataset == 'Robinson')

mapped_gray <<- greyfox %>% 
  filter(Type == 'Mapping') %>%
  filter(Dataset == 'Sacks')

paired <- greyfox %>% 
  filter(Type == 'ProperlyPaired') %>%
  filter(Dataset == 'Robinson')

paired_gray <- greyfox %>% 
  filter(Type == 'ProperlyPaired') %>%
  filter(Dataset == 'Sacks')

#Create summaries
paired_sum <- paired %>% group_by(Genome) %>%
  dplyr::summarize(Percent_Paired = mean(Score), sd_var1 = sd(Score, na.rm=TRUE))

paired_sum_gray <- paired_gray %>% group_by(Genome) %>%
  dplyr::summarize(Percent_Paired = mean(Score), sd_var1 = sd(Score, na.rm=TRUE))

mapped_sum <- mapped %>% group_by(Genome) %>%
  dplyr::summarize(Percent_Mapped = mean(Score), sd_var1 = sd(Score, na.rm=TRUE))

mapped_gray_sum <- mapped_gray %>% group_by(Genome) %>%
  dplyr::summarize(Percent_Mapped = mean(Score), sd_var1 = sd(Score, na.rm=TRUE))

depth_sum <- auto_depth %>% group_by(Genome) %>%
  dplyr::summarize(Avg_Depth = mean(Score), sd_var1 = sd(Score, na.rm=TRUE))

depth_gray_sum <- auto_depth_gray %>% group_by(Genome) %>%
  dplyr::summarize(Avg_Depth = mean(Score), sd_var1 = sd(Score, na.rm=TRUE))

#Convert from long to wide for T tests
#paired_wide <- spread(paired, Genome, Score)

test <- pairwise.wilcox.test(paired$Score, paired$Genome, p.adjust.method = "none")
test <- pairwise.wilcox.test(paired_gray$Score, paired_gray$Genome, p.adjust.method = "none")
test <- pairwise.wilcox.test(mapped$Score, mapped$Genome, p.adjust.method = "none")
test <- pairwise.wilcox.test(mapped_gray$Score, mapped_gray$Genome, p.adjust.method = "none")
test <- pairwise.wilcox.test(auto_depth$Score, auto_depth$Genome, p.adjust.method = "none")
test <- pairwise.wilcox.test(auto_depth_gray$Score, auto_depth_gray$Genome, p.adjust.method = "none")

#Plots
ggplot(data = depth, aes(Genome, as.numeric(Score), color = Dataset)) + 
  geom_violin() +
  theme_bw()

a<- ggplot(data = auto_depth, aes(Genome, as.numeric(Score), fill = Dataset)) + 
  geom_violin() +
  theme_bw() +
  theme(axis.text.x=element_text(color="black")) +
  xlab("") + ylab("Average autosome depth")

ggplot(data = mapped, aes(Genome, as.numeric(Score), color = Dataset)) + 
  geom_violin() +
  theme_bw()

b <- ggplot(data = paired, aes(Genome, as.numeric(Score), fill = Dataset)) + 
  geom_violin() +
  theme_bw() +
  theme(axis.text.x=element_text(color="black")) +
  xlab("") + ylab("Percent properly paired reads")

a/b
