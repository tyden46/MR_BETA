#Load libraries
library(plyr)
library(rentrez)
library(stringr)
library(tibble)
library(tidyselect)
library(ggplot2)
#Change working directory
setwd("C:/Users/tyden46/Desktop/MR_BETA_Tutorial")
#Read in the trueAbundance.txt file that has our organisms and their abundances
proportions=read.csv("trueAbundance.txt", sep="\t", header=FALSE)
listOfRefSeqIds=proportions$V1
listOfRefSeqIds=as.character(listOfRefSeqIds)
proportions=as.character(proportions$V2)
trueSpecies=c()
trueTaxonomy=c()
trueID=listOfRefSeqIds
fullEntryNameList=c()
genusBinom=c()
#Look up every ID in the NCBI nucleotide database to fetch their taxonomy
for (x in 1:length(trueID)){
  print(paste("Fetching species number ", x, "...", sep=" "))
  result=entrez_search(db="nuccore", term=trueID[x])
  textResult=entrez_fetch(db="nuccore", id=result$ids, rettype = "gb", retmode = "xml")
  myFile=XML::xmlToList(textResult)
  binomialNomenclature=myFile$GBSeq$GBSeq_organism
  taxonomy=myFile$GBSeq$GBSeq_taxonomy
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
  trueSpecies=append(trueSpecies, paste(myGenus, species, sep="-"))
  trueTaxonomy=append(trueTaxonomy, taxonomy)
}

#Create a list of true
taxOrder=c()
kingdom=c()
phylum=c()
class=c()
order=c()
family=c()
genus=c()
species=c()
q=1
#Organism lineages into arrays
for (x in trueTaxonomy){
  taxOrder=str_extract_all(x, "([A-Z]{1})(\\w+)")
  if (!is.null(taxOrder[[1]][1])){
    kingdom=append(kingdom, taxOrder[[1]][1])
  } else{
    kingdom=append(kingdom, "NA")
  }
  
  if (!is.null(taxOrder[[1]][2])){
    phylum=append(phylum, taxOrder[[1]][2])
  } else{
    phylum=append(phylum, "NA")
  }
  
  if (!is.null(taxOrder[[1]][3])){
    class=append(class, taxOrder[[1]][3])
  } else{
    class=append(class, "NA")
  }
  
  if (!is.null(taxOrder[[1]][4])){
    order=append(order, taxOrder[[1]][4])
  } else{
    order=append(order, "NA")
  }
  
  if (!is.null(taxOrder[[1]][5])){
    family=append(family, taxOrder[[1]][5])
  } else{
    family=append(family, "NA")
  }
  
  if (!is.null(taxOrder[[1]][6])){
    genus=append(genus, taxOrder[[1]][6])
  } else{
    genus=append(genus, "NA")
  }
  
  if (!is.null(taxOrder[[1]][7])){
    species=append(species, taxOrder[[1]][7])
  } else{
    species=append(species, "NA")
  }
  q=q+1
}
species=trueSpecies
#Create tibbles to hold lineages. Each row will be the lineage of an organism
taxClassSave=tibble(listOfRefSeqIds, kingdom, phylum, class, order, family, genus, species)
taxClassTrue=tibble(listOfRefSeqIds, kingdom, phylum, class, order, family, genus, species)
#Sometimes the genus isn't captured in the lineage properly so we look for NAs:
for (x in 1:length(taxClassTrue$listOfRefSeqIds)){
  if (is.na(taxClassTrue$genus[x])){
    taxClassTrue$genus[x]=genusBinom[[x]][1]
  }
}

#Manually review your lineages and try to replace any NA values with
#True taxa names
taxClassTrue$family[31]="Lachnospiraceae"
taxClassTrue$genus[31]="Anaerostipes"
taxClassTrue$species[31]="hadrus"

taxClassTrue$family[38]="unclassified Clostridiales"
taxClassTrue$genus[38]="unclassified Clostridiales (miscellaneous)"
taxClassTrue$species[38]="butyrate-producing bacterium"

taxClassTrue$family[40]="unclassified Clostridiales"
taxClassTrue$genus[40]="unclassified Clostridiales (miscellaneous)"
taxClassTrue$species[40]="butyrate-producing bacterium"

#Scale the abundance values so they sum to 1
taxClassTrue$proportion=proportions
taxClassTrue$proportion=(1/sum(as.numeric(taxClassTrue$proportion)))*as.numeric(taxClassTrue$proportion)

#Combine duplicates
taxClassTrue$combined <- apply( taxClassTrue[ , 2:8 ] , 1 , paste , collapse = "_" )
taxClassNoDuplicates=ddply(taxClassTrue, "combined", numcolwise(sum))
duplicates=duplicated(taxClassTrue[2:8])
idx=which(duplicates, arr.ind = TRUE)
taxClassTrue=taxClassTrue[-idx,]
final=merge(taxClassNoDuplicates, taxClassTrue, by="combined", all.x = FALSE)
final=final[2:10]

presentKingdoms=unique(taxClassTrue$kingdom)
presentPhyla=unique(taxClassTrue$phylum)
presentClass=unique(taxClassTrue$class)
presentOrder=unique(taxClassTrue$order)
presentFamily=unique(taxClassTrue$family)
presentGenus=unique(taxClassTrue$genus)
presentSpecies=unique(taxClassTrue$species)

#Create key value pairs for every taxa present and their abundance values
kingdomList=list()
for(x in presentKingdoms){
  key = x
  value = sum((final[which(final$kingdom==x), ])$proportion.x)
  kingdomList[[key]]=value
}

phylaList=list()
for(x in presentPhyla){
  key = x
  value = sum((final[which(final$phylum==x), ])$proportion.x)
  phylaList[[key]]=value
}

classList=list()
for(x in presentClass){
  key = x
  value = sum((final[which(final$class==x), ])$proportion.x)
  classList[[key]]=value
}

orderList=list()
for(x in presentOrder){
  key = x
  value = sum((final[which(final$order==x), ])$proportion.x)
  orderList[[key]]=value
}

familyList=list()
for(x in presentFamily){
  key = x
  value = sum((final[which(final$family==x), ])$proportion.x)
  familyList[[key]]=value
}

genusList=list()
for(x in presentGenus){
  key = x
  value = sum((final[which(final$genus==x), ])$proportion.x)
  genusList[[key]]=value
}

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

#Parse kraken's taxonomic classification report
krakenOutput=read.csv("krakenReport.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
krakenOutput2=krakenOutput
krakenOutput2$V1=(krakenOutput$V2/10000000)*100
krakenOutput=krakenOutput2
krakenOutput$V6=trimws(as.character(krakenOutput$V6), which="left")
#Need to change so that binomial nomenclature gets collapsed to just species name
speciesOnly=krakenOutput[which(krakenOutput$V4=="S"), ]
speciesNoGenus=c()
speciesFromKraken=speciesOnly$V6
for(x in 1:length(speciesFromKraken)){
  speciesName=paste(str_extract_all(speciesFromKraken[x], "[:alnum:]+")[[1]][1], str_extract_all(speciesFromKraken[x], "[:alnum:]+")[[1]][2], sep="-")
  speciesFromKraken[x]=speciesName
}
speciesOnly$V6=speciesFromKraken
krakenOutput[which(krakenOutput$V4=="S"), ]$V6=speciesFromKraken
krakenOutput=krakenOutput[-which(krakenOutput$V4=="S1"),]

krakenOutput$V1=as.numeric(krakenOutput$V1)/100

#Create key-value pairs for organims, their true abundance, and Kraken's estimate
krakenVsTruth=list()
for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="D"){
    if(!is.null(kingdomList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("Kingdom", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], kingdomList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("Kingdom", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}
for (x in names(kingdomList)){
  krakenKingdom=krakenOutput[which(krakenOutput$V4=="D"), ]
  idx=which(krakenKingdom$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("kingdom", x, sep="-")]]=c(0, kingdomList[[x]])
  }
}

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="P"){
    if(!is.null(phylaList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("phylum", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], phylaList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("phylum", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}

for (x in names(phylaList)){
  krakenPhylum=krakenOutput[which(krakenOutput$V4=="P"), ]
  idx=which(krakenPhylum$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("phylum", x, sep="-")]]=c(0, phylaList[[x]])
  }
}

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="C"){
    if(!is.null(classList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("class", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], classList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("class", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}
for (x in names(classList)){
  krakenClass=krakenOutput[which(krakenOutput$V4=="C"), ]
  idx=which(krakenClass$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("class", x, sep="-")]]=c(0, classList[[x]])
  }
}

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="O"){
    if(!is.null(orderList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("order", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], orderList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("order", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}
for (x in names(orderList)){
  krakenOrder=krakenOutput[which(krakenOutput$V4=="O"), ]
  idx=which(krakenOrder$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("order", x, sep="-")]]=c(0, orderList[[x]])
  }
}

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="F"){
    if(!is.null(familyList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("family", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], familyList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("family", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}

for (x in names(familyList)){
  krakenFamily=krakenOutput[which(krakenOutput$V4=="F"), ]
  idx=which(krakenFamily$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("family", x, sep="-")]]=c(0, familyList[[x]])
  }
}

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="G"){
    if(!is.null(genusList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("genus", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], genusList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("genus", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}
for (x in names(genusList)){
  krakenGenus=krakenOutput[which(krakenOutput$V4=="G"), ]
  idx=which(krakenGenus$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("genus", x, sep="-")]]=c(0, genusList[[x]])
  }
}

for (x in 1:length(krakenOutput$V1)){
  if(krakenOutput$V4[x]=="S"){
    if(!is.null(speciesList[[krakenOutput$V6[x]]])){
      krakenVsTruth[[paste("species", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], speciesList[[krakenOutput$V6[x]]])
    }
    else{
      krakenVsTruth[[paste("species", krakenOutput$V6[x], sep="-")]]=c(krakenOutput$V1[x], 0)
    }
  }
}
for (x in names(speciesList)){
  krakenSpecies=krakenOutput[which(krakenOutput$V4=="S"), ]
  idx=which(krakenSpecies$V6==x)
  if(identical(idx, integer(0))){
    krakenVsTruth[[paste("species", x, sep="-")]]=c(0, speciesList[[x]])
  }
}

#Coerce key-value pairs to a dataframe
dfTest <- ldply (krakenVsTruth, data.frame)
taxa=dfTest[seq(1, length(dfTest$.id),2 ),1 ]
krakenEstimate=dfTest[seq(1, length(dfTest$.id),2 ),2]
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
speciesIndex=c()
q=1
for (x in finalKrakenAndTruth$Taxa){
  if(str_detect(x, "species")){
    speciesIndex=append(speciesIndex, q)
  }
  q=q+1
}
finalKrakenAndTruthSave=finalKrakenAndTruth

#Identify overestimated and underestimated taxa
finalKrakenAndTruth=finalKrakenAndTruth[speciesIndex, ]
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

#Plot the overestimated species
g=ggplot(overestimated, aes(y=Ratio, x=reorder(Taxa, -Ratio))) +
  geom_col()+
  xlab("Taxa")+
  ylab(expression('Ratio'))+
  theme_grey()+
  labs(title="Overestimated Species",
       subtitle=paste("Ratio of Kraken2 Species Abundance Estimation", "to", "Actual Species Abundance", sep = "\n"),
       caption=str_wrap(paste("Fig 1.", "These species were correctly recognized by Kraken2
                           but their abundance values were overestimated. The figure shows the
                           species on the x-axis and the ratio of the Kraken2 estimate to actual
                           abundance on the y-axis.", sep = "\n")), width = 25)+
  theme(plot.title = element_text(hjust = 0.5, size=90),
        plot.subtitle = element_text(hjust = 0.5, size=40),
        axis.text = element_text(size=40, angle=90),
        axis.title = element_text(size=40), plot.caption = element_text(size=45, hjust=0))
ggsave("Overestimated.png", plot = last_plot(), device = "png", dpi = 400, height = 25, width = 25, limitsize = TRUE)

saveUnderestimated=underestimated
underestimated$Ratio=-log10(underestimated$Ratio)

#Plot underestimated species
g=ggplot(underestimated, aes(y=Ratio, x=reorder(Taxa, Ratio))) +
  geom_col()+
  xlab("Taxa")+
  ylab(expression('-log'[10]*'Ratio'))+
  theme_grey()+
  labs(title="Underestimated Species",
       subtitle=paste("Ratio of Actual Species Abundance", "to", "Kraken2 Species Abundance Estimation", sep = "\n"),
       caption=str_wrap(paste("Fig 2.", "These species were correctly recognized by Kraken2
                           but their abundance values were underestimated The figure shows the
                           species on the x-axis and the -log of the ratio of the Kraken2 estimate to actual
                           abundance on the y-axis.", sep = "\n")), width = 25)+
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
  labs(title="Misattributed Species",
       subtitle=paste("Relative Abundance of Species", "Identified in Kraken2 but not Present in Sample", sep = "\n"),
       caption=str_wrap(paste("Fig 3.", "These species were incorrectly recognized by Kraken2.
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
  labs(title="Unattributed Species",
       subtitle=paste("Relative Abundance of Species", "Not Identified in Kraken2 but Present in Sample", sep = "\n"),
       caption=str_wrap(paste("Fig 4.", "These species were not recognized by Kraken2
                           but were present in the sample. The figure shows the
                           species on the x-axis and the true relative abundance y-axis.", sep = "\n")), width = 25)+
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
score=sum(overestimated$Ratio)+sum(1/underestimated$Ratio)
write.csv(allOverEstimated[OESpeciesIndex, ], "overestimatedSpecies.csv")
write.csv(allUnderEstimated[UESpeciesIndex, ], "underestimatedSpecies.csv")

