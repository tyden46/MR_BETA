library(dplyr)
library(rentrez)
library(stringr)
library(tibble)
library(tidyselect)
library(ggplot2)

#Change Working Directory
setwd("C:/Users/tyden46/Desktop/MR_BETA_Tutorial")
#Read the previously generated "bias reports"
overEstimatedSpecies=read.csv("overestimatedSpecies.csv")
underEstimatedSpecies=read.csv("underestimatedSpecies.csv")
#Read the kraken report file containing Kraken 2's taxonomic classification
krakenOutput=read.csv("krakenReport.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
#Divide the second column by the total number of reads
krakenOutput$V1=krakenOutput$V2/10000000 #Change this depending on how many reads were generated
#Clean some of the spacing
krakenOutput$V6=trimws(as.character(krakenOutput$V6), which="left")
#Identify the entries which correspond to species and put them in a separate dataframe
speciesOnly=krakenOutput[which(krakenOutput$V4=="S"), ]
speciesOnly$V1=speciesOnly$V1*(1/(sum(speciesOnly$V1)))
speciesFromKraken=speciesOnly$V6
#Clean some of the species names
for(x in 1:length(speciesFromKraken)){
  speciesName=paste(str_extract_all(speciesFromKraken[x], "[:alnum:]+")[[1]][1], str_extract_all(speciesFromKraken[x], "[:alnum:]+")[[1]][2], sep="-")
  speciesFromKraken[x]=speciesName
}
#Substitute species names for cleaned species names
speciesOnly$V6=speciesFromKraken
krakenOutput[which(krakenOutput$V4=="S"), ]$V6=speciesFromKraken
#Scale so that species abundance sums to 1
speciesOnly$V1=(1/sum(speciesOnly$V1))*speciesOnly$V1
saveAb=speciesOnly$V1
#If any species are represented twice then find the sum of their abundances
combined=aggregate(speciesOnly$V1, by=list(V6=speciesOnly$V6), FUN=sum)
colnames(combined)=c("V6","V1")
#Eliminate duplicated species
speciesOnly=speciesOnly[-which(duplicated(speciesOnly$V6)),]
#Merge tables
speciesOnly=merge(combined,speciesOnly[, c(2:6)], by="V6", all.x = TRUE, all.y=FALSE)
speciesOnly$adjustedAbundance=speciesOnly$V1
#Identify unattributed species
for(x in underEstimatedSpecies$Taxa){
  query=str_remove(x, "species-")
  index=which(speciesOnly$V6==query)
  if(length(speciesOnly[index, ]$V6)<=0){
    de<-data.frame(query, 0, 0, 0, "S", 0, underEstimatedSpecies[which(underEstimatedSpecies$Taxa==x), ]$True.Abundance)
    names(de)<-c("V6","V1", "V2", "V3", "V4", "V5", "adjustedAbundance")
    speciesOnly <- rbind(speciesOnly, de)
  }
}
#For all of the species in Kraken 2's estimate, find which are overestimated, underestimated
#Or Misattributed and scale down, scale up, or remove accordingly
for (x in speciesOnly$V6){
  query=paste("species", x, sep="-")
  index=which(overEstimatedSpecies$Taxa==query)
  index2=which(underEstimatedSpecies$Taxa==query)
  #If the species is in our bias profile for overestimated species
  if(length(overEstimatedSpecies[index, ]$True.Abundance)>0){
    #If the species is misattributed (ie shouldn't be in the sample)
    if(overEstimatedSpecies[index, ]$True.Abundance==0){
      speciesIndex=which(speciesOnly$V6==x)
      speciesOnly[which(speciesOnly$V6==x), ]$adjustedAbundance=0
      saveAb=saveAb[-speciesIndex]
    } else if(overEstimatedSpecies[index, ]$True.Abundance>0){ #If the Species is present but overestimated
      speciesIndex=which(speciesOnly$V6==x)
      speciesOnly[which(speciesOnly$V6==x), ]$adjustedAbundance=(speciesOnly[which(speciesOnly$V6==x), ]$adjustedAbundance)/(overEstimatedSpecies[index, ]$Ratio)
      saveAb=saveAb[-speciesIndex]
    }
  }
  #If the species is in our bias profile for underestimated species
  if(length(underEstimatedSpecies[index2, ]$True.Abundance)>0){
    if(underEstimatedSpecies[index2, ]$Ratio>0){
      speciesIndex=which(speciesOnly$V6==x)
      speciesOnly[which(speciesOnly$V6==x), ]$adjustedAbundance=(speciesOnly[which(speciesOnly$V6==x), ]$adjustedAbundance)/(underEstimatedSpecies[index2, ]$Ratio)
      saveAb=saveAb[-speciesIndex]
    }
  }
  #If the species is misattributed
  if(length(underEstimatedSpecies[index, ]$True.Abundance)<=0 && length(underEstimatedSpecies[index2, ]$True.Abundance)<=0){
    speciesOnly[which(speciesOnly$V6==x), ]$adjustedAbundance=0
  }
}

#Scale so that the abundance values sum to one
speciesOnly$adjustedAbundance=(1/(sum(speciesOnly$adjustedAbundance)))*speciesOnly$adjustedAbundance
merged=rbind(overEstimatedSpecies, underEstimatedSpecies)
merged$Taxa=as.character(merged$Taxa)
#Clean some of the species names
for (x in 1:length(merged$Taxa)){
  merged$Taxa[x]=as.character(str_remove(merged$Taxa[x], "species-"))
}
#Clean and join tables
colnames(speciesOnly)[1]="Taxa"
speciesOnly=speciesOnly[,c(1,2,7)]
merged=merged[order(merged$Taxa), ]
final=right_join(merged, speciesOnly)#, by="Taxa", all.x = TRUE, all.y = TRUE)
final[is.na(final)] <- 0
#Scale so that adjusted abundance sums to 1
final$adjustedAbundance=(1/(sum(final$adjustedAbundance)))*final$adjustedAbundance
#Scale so that true abundance sums to 1
final$True.Abundance=(1/(sum(final$True.Abundance)))*final$True.Abundance
final$postDifference=final$adjustedAbundance - final$True.Abundance
final$difference=final$V1 - final$True.Abundance
before=c("Before correction", sum(abs(final$difference)))

#Graph results
finalResult=c("After correction", sum(abs(final$postDifference)))
beforeAndAfter=as.data.frame(rbind(before,finalResult), stringsAsFactors = FALSE)
beforeAndAfter$V1=factor(beforeAndAfter$V1, levels = c("Before correction", "After correction"))
g=ggplot(beforeAndAfter, aes(y=as.numeric(V2), x=V1)) +
  geom_col()+
  xlab("Step in Correction Process")+
  ylab("Difference between True and Estimated Abundance")+
  theme_grey()+
  labs(title="Reassignment Results",
       subtitle=paste("Bias adjustment algorithm increases", "accuracy after each step", sep = "\n"),
       caption=str_wrap(paste("Fig 5.", "Taxonomic profile adjustment based on the bias report
                              increases accuracy In total, error rate was reduced from", sum(abs(final$difference)), " to ", sum(abs(final$postDifference)), sep = "\n")), width = 25)+
  theme(plot.title = element_text(hjust = 0.5, size=90),
        plot.subtitle = element_text(hjust = 0.5, size=40),
        axis.text = element_text(size=40, angle=90),
        axis.title = element_text(size=40), plot.caption = element_text(size=45, hjust=0))
ggsave("correction.png", plot = last_plot(), device = "png", dpi = 400, height = 25, width = 25, limitsize = TRUE)
final$diffdiff=abs(final$difference-final$postDifference)
present=final[which(final$adjustedAbundance>0), ]
output=present[,c(2,8)]
write.csv(output, "AdjustedOutput.csv", row.names = FALSE)