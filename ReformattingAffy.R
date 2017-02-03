#Reformatting the Affy data

#setwd to Affy

library(plyr)
library(affy)
library(multtest)

#Affy files are located:
setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")

#For testing:
BrainRegion<-"ACG"
BrainRegionFileName<-"ACG"
AnnotationColumns<-
  
############CODE for function to reformat LCM_AMY results into generalized output######################

Reformat_Affy<-function(BrainRegion, BrainRegionFileName){
  
  LM4Betas<-read.table(paste(BrainRegionFileName, "_ModelBetas.csv", sep=""), header=T, row.names=1, sep=",")
  #Hmm... note: I'm not going to be able to automate all of this easily, because only some regions have "EarlyLateCohort" as a variable - although those results may not be worth extracting anyway, right?  Intercept probably isn't worth extracting either...I don't think many of these analyses had pH or Age centered when they were run.
  
  colnames(LM4Betas)
  head(LM4Betas)
  str(LM4Betas)
  
  LM4BetasExpression<-LM4Betas[,c(4:11)]
  colnames(LM4BetasExpression)
  rm(LM4Betas)
  
  #I'm a little nervous about automatically renaming these - I should probably print something to the screen so that I double check that it is done properly.
  
  VariablesOfInterest<-c("BP", "MDD", "SCHIZ", "BrainPH", "AgonalFactor", "PMI",  "Gender", "Age")
  BetaColumnNames<-VariablesOfInterest
  
  for(i in 1:8){
    BetaColumnNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "Beta", sep="_")
  }
  
  print(cbind(colnames(LM4BetasExpression), BetaColumnNames))
  
  colnames(LM4BetasExpression)<-BetaColumnNames
  
  
  LM4DirectionOfEffect<-matrix("Up", nrow(LM4BetasExpression), ncol(LM4BetasExpression))
  
  LM4DirectionOfEffect<-LM4BetasExpression
  LM4DirectionOfEffect[LM4DirectionOfEffect>0]<-"Up"
  LM4DirectionOfEffect[LM4DirectionOfEffect<0]<-"Down"
  #table(LM4DirectionOfEffect[,3])
  #head(cbind(LM4DirectionOfEffect,LM4BetasExpression))
  
  colnames(LM4DirectionOfEffect)
  for(i in 1:8){
    colnames(LM4DirectionOfEffect)[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "DirectionOfEffect", sep="_")
  }
  
  
  LM4Tstat<-read.table(paste(BrainRegionFileName, "_ModelTstat.csv", sep=""), header=T, row.names=1, sep=",")
  colnames(LM4Tstat)
  LM4TstatExpression<-LM4Tstat[,c(4:11)]
  rm(LM4Tstat)
  
  TstatColumnNames<-VariablesOfInterest
  for(i in 1:8){
    TstatColumnNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "Tstat", sep="_")
  }
  
  print(cbind(colnames(LM4TstatExpression), TstatColumnNames))
  
  colnames(LM4TstatExpression)<-TstatColumnNames
  
  LM4pvalues<-read.table(paste(BrainRegionFileName, "_ModelpvaluesRAW.csv", sep=""), header=T, row.names=1, sep=",")
  colnames(LM4pvalues)
  ProbeInfo<-LM4pvalues[,c(1:2)]
  
  LM4pvaluesExpression<-LM4pvalues[,c(4:11)]
  rm(LM4pvalues)
  
  PvalueColumnNames<-VariablesOfInterest
  for(i in 1:8){
    PvalueColumnNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "NominalPvalue", sep="_")
  }
  
  print(cbind(colnames(LM4pvaluesExpression), PvalueColumnNames))
  
  colnames(LM4pvaluesExpression)<-PvalueColumnNames
  
  PercentileRank<-LM4pvaluesExpression
  
  for(i in c(1:ncol(LM4pvaluesExpression))){
    PercentileRank[,i]<-rank(as.numeric(LM4pvaluesExpression[,i]), ties.method=("average"))/length(LM4pvaluesExpression[,i])
  }
  
  for(i in 1:8){
    colnames(PercentileRank)[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "PercentileRank", sep="_")
  }
  
  head(cbind(LM4pvaluesExpression, PercentileRank))
  
  
  #All of this needs to be changed - the adjusted pvalues for diagnosis can be read in, but all of the rest will need to be recalculated.

  
  #I need to make sure I'm following the same order as VariablesOfInterest
 
  colnames(ProbeInfo)<-c("GeneSetUnigeneClusterID", "SYMBOLREANNOTATED")
  ProbeAnnotation<-ProbeInfo
  colnames(ProbeAnnotation)
  head(ProbeAnnotation)
  
  
  ProbeAnnotationColNames<-colnames(ProbeAnnotation)
  
  for(i in 1:ncol(ProbeAnnotation)){
    ProbeAnnotationColNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ProbeAnnotation", colnames(ProbeAnnotation)[i], sep="_")
  }
  
  print(cbind(ProbeAnnotationColNames, colnames(ProbeAnnotation)))
  
  colnames(ProbeAnnotation)<-ProbeAnnotationColNames
  
  BH_AdjustedPval_FDR<-LM4pvaluesExpression
  
  #Applying multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
  
  for(i in c(1:8)){
  TempPvalAdj<-mt.rawp2adjp(LM4pvaluesExpression[,i], proc="BH")
  PvalAdj<-TempPvalAdj$adjp[order(TempPvalAdj$index),]
  BH_AdjustedPval_FDR[,i]<-PvalAdj[,2]
  }
  
  head(BH_AdjustedPval_FDR)
  
  for(i in 1:8){
    colnames(BH_AdjustedPval_FDR)[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "FDR", sep="_")
  }
  
  temp<<-data.frame(ProbeAnnotation,LM4BetasExpression, LM4DirectionOfEffect, LM4TstatExpression, LM4pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
  
}

############End of CODE for function to reformat Affy results into generalized output######################

rm(BH_AdjustedPval_FDR, Affy_Ref8v2Annotation_fromGEO, LM4BetasExpression, LM4DirectionOfEffect, LM4pvaluesExpression, LM4TstatExpression, PercentileRank, ProbeInfo, ProbeAnnotation, PvalAdj, temp)
rm(BetaColumnNames, ProbeAnnotationColNames, TempPvalAdj, TstatColumnNames, PvalueColumnNames, VariablesOfInterest)


#Change to output working directory

setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")
Reformat_Affy("ACG", "ACg")
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "ACG_Affymetrix_Affy6Region_ModelLM4.csv")

setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")
Reformat_Affy("AMY", "Amy")
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "AMY_Affymetrix_Affy6Region_ModelLM4.csv")

setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")
Reformat_Affy("CB", "CB")
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "CB_Affymetrix_Affy6Region_ModelLM4.csv")

setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")
Reformat_Affy("DLPFC", "DLPFC")
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "DLPFC_Affymetrix_Affy6Region_ModelLM4.csv")

setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")
Reformat_Affy("HC", "HC")
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "HC_Affymetrix_Affy6Region_ModelLM4.csv")

setwd("~/Documents/Affy/ElyseReanalysis/Affy Summary for Database")
Reformat_Affy("NACC", "NAcc")
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "NACC_Affymetrix_Affy6Region_ModelLM4.csv")


##################
#The DLPFC required special treatment because it is my output, not Elyse's standardized output.

#For testing:
BrainRegion<-"DLPFC"
BrainRegionFileName<-"DLPFC"
AnnotationColumns<-
  
  setwd("~/Documents/Affy/NoPC1correct/DLPFC circadian/LM4_basic")

LM4Betas<-read.table("GeneByCellTypeSubjVar2_Betas.csv", header=T, row.names=1, sep=",")
#Hmm... note: I'm not going to be able to automate all of this easily, because only some regions have "EarlyLateCohort" as a variable - although those results may not be worth extracting anyway, right?  Intercept probably isn't worth extracting either...I don't think many of these analyses had pH or Age centered when they were run.

colnames(LM4Betas)
head(LM4Betas)
str(LM4Betas)

LM4BetasExpression<-LM4Betas[,c(2:9)]
colnames(LM4BetasExpression)
rm(LM4Betas)

#I'm a little nervous about automatically renaming these - I should probably print something to the screen so that I double check that it is done properly.

VariablesOfInterest<-c("BrainPH", "AgonalFactor", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ")
BetaColumnNames<-VariablesOfInterest

for(i in 1:8){
  BetaColumnNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "Beta", sep="_")
}

print(cbind(colnames(LM4BetasExpression), BetaColumnNames))

colnames(LM4BetasExpression)<-BetaColumnNames


LM4DirectionOfEffect<-matrix("Up", nrow(LM4BetasExpression), ncol(LM4BetasExpression))

LM4DirectionOfEffect<-LM4BetasExpression
LM4DirectionOfEffect[LM4DirectionOfEffect>0]<-"Up"
LM4DirectionOfEffect[LM4DirectionOfEffect<0]<-"Down"
table(LM4DirectionOfEffect[,3])
head(cbind(LM4DirectionOfEffect,LM4BetasExpression))

colnames(LM4DirectionOfEffect)
for(i in 1:8){
  colnames(LM4DirectionOfEffect)[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "DirectionOfEffect", sep="_")
}


LM4Tstat<-read.table("GeneByCellTypeSubjVar2_Tstat.csv", header=T, row.names=1, sep=",")
colnames(LM4Tstat)
LM4TstatExpression<-LM4Tstat[,c(2:9)]
rm(LM4Tstat)

TstatColumnNames<-VariablesOfInterest
for(i in 1:8){
  TstatColumnNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "Tstat", sep="_")
}

print(cbind(colnames(LM4TstatExpression), TstatColumnNames))

colnames(LM4TstatExpression)<-TstatColumnNames

LM4pvalues<-read.table("GeneByCellTypeSubjVar2_Pvalues.csv", header=T, row.names=1, sep=",")
colnames(LM4pvalues)
ProbeInfo<-LM4pvalues[,c(10:11)]

LM4pvaluesExpression<-LM4pvalues[,c(2:9)]
rm(LM4pvalues)

PvalueColumnNames<-VariablesOfInterest
for(i in 1:8){
  PvalueColumnNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "NominalPvalue", sep="_")
}

print(cbind(colnames(LM4pvaluesExpression), PvalueColumnNames))

colnames(LM4pvaluesExpression)<-PvalueColumnNames

PercentileRank<-LM4pvaluesExpression

for(i in c(1:ncol(LM4pvaluesExpression))){
  PercentileRank[,i]<-rank(as.numeric(LM4pvaluesExpression[,i]), ties.method=("average"))/length(LM4pvaluesExpression[,i])
}

for(i in 1:8){
  colnames(PercentileRank)[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "PercentileRank", sep="_")
}

head(cbind(LM4pvaluesExpression, PercentileRank))


#All of this needs to be changed - the adjusted pvalues for diagnosis can be read in, but all of the rest will need to be recalculated.


#I need to make sure I'm following the same order as VariablesOfInterest

colnames(ProbeInfo)<-c("GeneSetUnigeneClusterID", "SYMBOLREANNOTATED")
ProbeAnnotation<-ProbeInfo
colnames(ProbeAnnotation)
head(ProbeAnnotation)


ProbeAnnotationColNames<-colnames(ProbeAnnotation)

for(i in 1:ncol(ProbeAnnotation)){
  ProbeAnnotationColNames[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ProbeAnnotation", colnames(ProbeAnnotation)[i], sep="_")
}

print(cbind(ProbeAnnotationColNames, colnames(ProbeAnnotation)))

colnames(ProbeAnnotation)<-ProbeAnnotationColNames

BH_AdjustedPval_FDR<-LM4pvaluesExpression

#Applying multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):

for(i in c(1:8)){
  TempPvalAdj<-mt.rawp2adjp(LM4pvaluesExpression[,i], proc="BH")
  PvalAdj<-TempPvalAdj$adjp[order(TempPvalAdj$index),]
  BH_AdjustedPval_FDR[,i]<-PvalAdj[,2]
}

head(BH_AdjustedPval_FDR)

for(i in 1:8){
  colnames(BH_AdjustedPval_FDR)[i]<-paste(BrainRegion, "Affymetrix_Affy6Region_ModelLM4", VariablesOfInterest[i], "FDR", sep="_")
}

temp<<-data.frame(ProbeAnnotation,LM4BetasExpression, LM4DirectionOfEffect, LM4TstatExpression, LM4pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)



setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Affymetrix_Affy6Region")
write.csv(temp, "DLPFC_Affymetrix_Affy6Region_ModelLM4.csv")

