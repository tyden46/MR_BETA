#Load libraries
library(plyr)
library(rentrez)
library(stringr)
library(tibble)
library(tidyselect)
library(ggplot2)
library(scales)
#Change working directory
set_entrez_key("c20bee62610f497c2766ffd7e8430abc0408")
setwd("C:/Users/tyden46/Desktop/MR_BETA_Tutorial")
#Read in the trueAbundance.txt file that has our organisms and their abundances
proportions=read.csv("trueAbundanceTenThousandReads.txt", sep="\t", header=FALSE)
listOfRefSeqIds=proportions$V1
listOfRefSeqIds=as.character(listOfRefSeqIds)
proportions=as.character(proportions$V2)
trueSpecies=c()
trueTaxonomy=c()
trueID=listOfRefSeqIds
fullEntryNameList=c()
genusBinom=c()
genomeSizeList=c()
#Look up every ID in the NCBI nucleotide database to fetch their taxonomy
for (x in 1:length(trueID)){
  print(paste("Fetching species number ", x, "...", sep=" "))
  result=entrez_search(db="nuccore", term=trueID[x])
  textResult=entrez_fetch(db="nuccore", id=result$ids, rettype = "gb", retmode = "xml")
  myFile=XML::xmlToList(textResult)
  binomialNomenclature=myFile$GBSeq$GBSeq_organism
  taxonomy=myFile$GBSeq$GBSeq_taxonomy
  genomeSize=as.numeric(myFile$GBSeq$GBSeq_length)
  genomeSizeList[[trueID[x]]]=genomeSize
  #Regular expressions can be used to fetch the lineage
  fullEntryName=str_extract_all(binomialNomenclature, "\\w+")
  genus=fullEntryName[1]
  genusBinom=append(genusBinom, genus)
  fullEntryNameList=append(fullEntryNameList, fullEntryName)
  binomialNomRemoveSubsp=str_extract(binomialNomenclature, "(([A-z]+)([:space:]{1})([A-z]+))")
  binomialNomRemoveSubsp=as.character(binomialNomRemoveSubsp)
  myGenus= str_extract(binomialNomRemoveSubsp, "([A-z]+){1}")
  myGenus=str_remove_all(myGenus, "\\[")
  myGenus=str_remove_all(myGenus, "\\]")
  species = str_extract(binomialNomRemoveSubsp, "[:space:]{1}([A-z]+){1}")
  species = str_replace(species, " ", "")
  trueSpecies=append(trueSpecies, paste(myGenus, species, sep=" "))
  trueTaxonomy=append(trueTaxonomy, taxonomy)
}
#Clean Up
rm(myFile)
rm(result)
rm(fullEntryName)
rm(genus)
rm(binomialNomenclature)
rm(binomialNomRemoveSubsp)
rm(genomeSize)
rm(myGenus)
rm(species)
rm(taxonomy)
rm(textResult)
rm(x)

#Create a list of true genuses and species
genus=c()
species=c()
q=1
#Organism lineages into arrays
for (x in trueTaxonomy){
  taxOrder=str_extract_all(x, "([A-Z]{1})(\\w+)")
  if (!is.null(taxOrder[[1]][6])){
    genus=append(genus, taxOrder[[1]][6])
  } else{
    genus=append(genus, "NA")
  }
  q=q+1
}
#Clean up
rm(taxOrder)
rm(x)
species=trueSpecies
#Create tibbles to hold lineages. Each row will be the lineage of an organism
taxClassSave=tibble(listOfRefSeqIds, genus, species)
taxClassTrue=tibble(listOfRefSeqIds, genus, species)
genomeSize=c(rep(0, length(taxClassSave$listOfRefSeqIds)))
ind=1
for (x in taxClassSave$listOfRefSeqIds){
  genomeSize[ind]=genomeSizeList[x]
  ind=ind+1
}
taxClassTrue$genomeSize=genomeSize
#Sometimes the genus isn't captured in the lineage properly so we look for NAs:
for (x in 1:length(taxClassTrue$listOfRefSeqIds)){
  if (is.na(taxClassTrue$genus[x])){
    taxClassTrue$genus[x]=genusBinom[[x]][1]
  }
}

#Manually review your lineages and try to replace any NA values with
#True taxa names
taxClassTrue$species[31]="butyrate-producing bacterium"

taxClassTrue$species[38]="butyrate-producing bacterium"

taxClassTrue$species[40]="butyrate-producing bacterium"

#Clean up
rm(ind)
rm(x)
rm(genomeSize)
rm(genomeSizeList)
rm(genus)
rm(listOfRefSeqIds)
rm(q)
rm(species)
rm(trueID)
rm(trueSpecies)
rm(trueTaxonomy)

#Scale the abundance values so they sum to 1
taxClassTrue$proportion=proportions
taxClassBeforeNormalization=taxClassTrue
newProportion=(as.numeric(taxClassTrue$proportion))/taxClassTrue$genomeSize
taxClassTrue$proportion=(1/sum(as.numeric(newProportion)))*as.numeric(newProportion)
#Combine duplicates
taxClassNoDuplicates=ddply(taxClassTrue, "species", numcolwise(sum))
duplicates=duplicated(taxClassTrue[3])
idx=which(duplicates, arr.ind = TRUE)
taxClassTrue=taxClassTrue[-idx,]
final=merge(taxClassNoDuplicates, taxClassTrue, by="species", all.x = FALSE)
final=final[c(1,3,4)]

presentSpecies=taxClassTrue$species

#Create key value pairs for every taxa present and their abundance values

speciesList=list()
for(x in presentSpecies){
  key = x
  value = sum((final[which(final$species==x), ])$proportion.x)
  speciesList[[key]]=value
}
sum=0
for (x in speciesList){
  sum=sum+x
}
rm(sum)

#Parse kraken's taxonomic classification report
krakenOutput=read.csv("C://Users//tyden46//Desktop//krakenReport-TenThousandReads-1.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
krakenOutput2=krakenOutput
krakenOutput2$V1=(krakenOutput$V2/10000)*100
krakenOutput=krakenOutput2
krakenOutput$V6=trimws(as.character(krakenOutput$V6), which="left")
#Need to change so that binomial nomenclature gets collapsed to just species name
speciesOnly=krakenOutput[which(krakenOutput$V4=="S"), ]
speciesNoGenus=c()
speciesFromKraken=speciesOnly$V6
for(x in 1:length(speciesFromKraken)){
  speciesName=paste(str_extract_all(speciesFromKraken[x], "[:alnum:]+")[[1]][1], str_extract_all(speciesFromKraken[x], "[:alnum:]+")[[1]][2], sep=" ")
  speciesFromKraken[x]=speciesName
}
speciesOnly$V6=speciesFromKraken
krakenOutput[which(krakenOutput$V4=="S"), ]$V6=speciesFromKraken
krakenOutput=krakenOutput[-which(krakenOutput$V4=="S1"),]

krakenOutput$V1=as.numeric(krakenOutput$V1)/100

#Create key-value pairs for organims, their true abundance, and Kraken's estimate
krakenVsTruth=list()

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="S"){
    if(!is.null(speciesList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[krakenOutput$V6[x]]]=c(krakenOutput$V1[x], speciesList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[krakenOutput$V6[x]]]=c(krakenOutput$V1[x], 0)
    }
  }
}
for (x in names(speciesList)){
  krakenSpecies=krakenOutput[which(krakenOutput$V4=="S"), ]
  idx=which(krakenSpecies$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[x]]=c(0, speciesList[[x]])
  }
}

#Coerce key-value pairs to a dataframe
dfTest <- ldply (krakenVsTruth, data.frame)
taxa=dfTest[seq(1, length(dfTest$.id),2 ),1 ]
krakenEstimate=dfTest[seq(1, length(dfTest$.id),2 ),2]
krakenEstimate=krakenEstimate*(1/sum(krakenEstimate))
trueAbundance=dfTest[seq(2, length(dfTest$.id),2 ),2]
#Calculate ratio between true and estimated abundance
ratio=krakenEstimate/trueAbundance

#Make dataframe pretty
finalKrakenAndTruth=tibble(taxa, as.numeric(krakenEstimate), as.numeric(trueAbundance), as.numeric(ratio))
finalKrakenAndTruth=dplyr::rename(finalKrakenAndTruth, Taxa = taxa)
finalKrakenAndTruth=dplyr::rename(finalKrakenAndTruth, `True Abundance`=`as.numeric(trueAbundance)`)
finalKrakenAndTruth=dplyr::rename(finalKrakenAndTruth, `Kraken Estimated Abundance`=`as.numeric(krakenEstimate)`)
finalKrakenAndTruth=dplyr::rename(finalKrakenAndTruth, `Ratio`=`as.numeric(ratio)`)
#Figure out where the species are in the dataframe
#We will only save the overestimated and underestimated species, not other taxa
finalKrakenAndTruthSave=finalKrakenAndTruth

#Identify overestimated and underestimated taxa
overestimated=finalKrakenAndTruth[which(finalKrakenAndTruth$Ratio>1), ]
inKrakenNotSample=finalKrakenAndTruth[which(is.infinite(finalKrakenAndTruth$Ratio)),]
overestimated=overestimated[which(is.finite(overestimated$Ratio)),]

underestimated=finalKrakenAndTruth[which(finalKrakenAndTruth$Ratio<1), ]
inSampleNotKraken=underestimated[which(underestimated$Ratio==0), ]
underestimated=underestimated[which(underestimated$Ratio>0), ]
underestimated=underestimated[which(is.finite(underestimated$Ratio)),]
underestimated$Ratio=1/underestimated$Ratio
#Doublecheck to remove any subspecies
idx=which(underestimated$Taxa=="species-sp")
if(length(idx)>0){
  underestimated=underestimated[-which(underestimated$Taxa=="species-sp"),]
}


overestimated=overestimated[order(overestimated$Ratio), ]
underestimated=underestimated[order(underestimated$Ratio), ]

j=1
myBreaks=c()
while(j<max(overestimated$Ratio+0.6)){
  myBreaks=append(myBreaks, j)
  j=j+0.4
}

#Plot the overestimated species
g=ggplot(overestimated, aes(y=Ratio, x=reorder(Taxa, -Ratio))) +
  geom_col()+
  xlab("Taxa")+
  scale_y_continuous(breaks=myBreaks, limits=c(1,max(overestimated$Ratio)+0.1),oob = rescale_none)+
  ylab(expression('Ratio'))+
  theme_grey()+
  labs(title="Kraken2 Overestimated Species",
       subtitle=paste("Ratio of Kraken2 Species Abundance Estimation", "to", "Actual Species Abundance", sep = "\n"),
       caption=str_wrap(paste("Fig 2.", "These species were correctly recognized by Kraken2
                           but their abundance values were overestimated. The figure shows the
                           species on the x-axis and the ratio of the Kraken2 estimate to actual
                           abundance on the y-axis.", "The total number of overestimated species is ", length(overestimated$Taxa), sep = "\n")), width = 25)+
  theme(plot.title = element_text(hjust = 0.5, size=90),
        plot.subtitle = element_text(hjust = 0.5, size=40),
        axis.text = element_text(size=40, angle=90),
        axis.title = element_text(size=40), plot.caption = element_text(size=45, hjust=0))
ggsave("Overestimated.png", plot = last_plot(), device = "png", dpi = 400, height = 25, width = 25, limitsize = TRUE)

saveUnderestimated=underestimated
k=-1
myBreaks2=c()
while(k>-max(underestimated$Ratio)){
  myBreaks2=append(myBreaks2, k)
  k=k-0.4
}

#Plot underestimated species
g=ggplot(underestimated, aes(y=-Ratio, x=reorder(Taxa, -Ratio))) +
  geom_col()+
  xlab("Taxa")+
  ylab(expression('Ratio'))+
  scale_y_continuous(breaks=myBreaks2, labels=(myBreaks2*-1), limits=c(-max(underestimated$Ratio)-0.1, -1),oob = rescale_none)+
  theme_grey()+
  labs(title="Kraken2 Underestimated Species",
       subtitle=paste("Ratio of Actual Species Abundance", "to", "Kraken2 Species Abundance Estimation", sep = "\n"),
       caption=str_wrap(paste("Fig 4.", "These species were correctly recognized by Kraken2
                           but their abundance values were underestimated The figure shows the
                           species on the x-axis and the -log of the ratio of the Kraken2 estimate to actual
                           abundance on the y-axis.", "The total number of underestimated species is ", length(underestimated$Taxa), sep = "\n")), width = 25)+
  theme(plot.title = element_text(hjust = 0.5, size=90),
        plot.subtitle = element_text(hjust = 0.5, size=40),
        axis.text = element_text(size=40, angle=90),
        axis.title = element_text(size=40), plot.caption = element_text(size=45, hjust=0))
ggsave("Underestimated.png", plot = last_plot(), device = "png", dpi = 400, height = 25, width = 25, limitsize = TRUE)

#Plot Misattributed Species
g=ggplot(inKrakenNotSample[order(inKrakenNotSample$`Kraken Estimated Abundance`, decreasing = TRUE)[1:20], ], aes(y=`Kraken Estimated Abundance`, x=reorder(Taxa, -`Kraken Estimated Abundance`))) +
  geom_col()+
  xlab("Taxa")+
  ylab(expression('Relative Abundance'))+
  theme_grey()+
  labs(title="Kraken2 Misattributed Species",
       subtitle=paste("Relative Abundance of Species", "Identified in Kraken2 but not Present in Sample", sep = "\n"),
       caption=str_wrap(paste("Fig 6.", "These species were incorrectly recognized by Kraken2.
       In truth, the abundance value of each should be 0. The figure shows the
       species on the x-axis and the abundance estimate of the Kraken2 on the y-axis. Only the top 20 species, sorted by estimated abundance, are shown. In total the number of misattributed species is: ", length(inKrakenNotSample$Taxa), sep = "\n")), width = 25)+
  theme(plot.title = element_text(hjust = 0.5, size=90),
        plot.subtitle = element_text(hjust = 0.5, size=40),
        axis.text = element_text(size=40, angle=90),
        axis.title = element_text(size=40), plot.caption = element_text(size=45, hjust=0))
ggsave("Misattributed.png", plot = last_plot(), device = "png", dpi = 400, height = 25, width = 25, limitsize = TRUE)

g=ggplot(inSampleNotKraken, aes(y=`True Abundance`, x=reorder(Taxa, -`True Abundance`))) +
  geom_col()+
  xlab("Taxa")+
  ylab(expression('Relative Abundance'))+
  theme_grey()+
  labs(title="Kraken2 Unattributed Species",
       subtitle=paste("Relative Abundance of Species", "Not Identified in Kraken2 but Present in Sample", sep = "\n"),
       caption=str_wrap(paste("Fig 8.", "These species were not recognized by Kraken2
                           but were present in the sample. The figure shows the
                           species on the x-axis and the true relative abundance y-axis.", "The total number of misattributed species is ", length(inSampleNotKraken$Taxa), sep = "\n")), width = 25)+
  theme(plot.title = element_text(hjust = 0.5, size=90),
        plot.subtitle = element_text(hjust = 0.5, size=40),
        axis.text = element_text(size=40, angle=90),
        axis.title = element_text(size=40), plot.caption = element_text(size=45, hjust=0))
ggsave("Unattributed.png", plot = last_plot(), device = "png", dpi = 400, height = 25, width = 25, limitsize = TRUE)

#Identify over- and under- estimated species and write to csv
#This constitutes the "bias profile" generated by MR BETA
allOverEstimated=finalKrakenAndTruth[which(finalKrakenAndTruth$Ratio>1), ]
allUnderEstimated=finalKrakenAndTruth[which(finalKrakenAndTruth$Ratio<1), ]
OESpeciesIndex=c()
for (x in 1:length(allOverEstimated$Taxa)){
  if(str_detect(allOverEstimated$Taxa[x], "species")){
    OESpeciesIndex=append(OESpeciesIndex, x)
  }
}

UESpeciesIndex=c()
for (x in 1:length(allUnderEstimated$Taxa)){
  if(str_detect(allUnderEstimated$Taxa[x], "species")){
    UESpeciesIndex=append(UESpeciesIndex, x)
  }
}
allOverEstimated$difference=allOverEstimated$`Kraken Estimated Abundance` - allOverEstimated$`True Abundance`
allUnderEstimated$difference=allUnderEstimated$`Kraken Estimated Abundance` - allUnderEstimated$`True Abundance`
write.csv(allOverEstimated[OESpeciesIndex, ], "overestimatedSpecies.csv")
write.csv(allUnderEstimated[UESpeciesIndex, ], "underestimatedSpecies.csv")

