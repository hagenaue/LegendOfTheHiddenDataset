
#setwd to LM2_OutputBest
#list.files()

#test run:
BrainRegion<-"AAA"


############CODE for function to reformat LCM_AMY results into generalized output######################

ReformatLCM_AMY<-function(BrainRegion, AnnotationColumns){
  
  LM2Betas<-read.table("LM2Betas.csv", header=T, row.names=1, sep=",")
  #Hmm... note: I'm not going to be able to automate all of this easily, because only some regions have "EarlyLateCohort" as a variable - although those results may not be worth extracting anyway, right?  Intercept probably isn't worth extracting either...I don't think many of these analyses had pH or Age centered when they were run.
  
  colnames(LM2Betas)
  head(LM2Betas)
  str(LM2Betas)
  
  LM2BetasExpression<-LM2Betas[,c(2:7)]
  colnames(LM2BetasExpression)
  rm(LM2Betas)
  
  #I'm a little nervous about automatically renaming these - I should probably print something to the screen so that I double check that it is done properly.
  
  VariablesOfInterest<-c("MDD", "Age", "RNAIntegrity", "BrainPH", "Gender", "PMI")
  BetaColumnNames<-VariablesOfInterest
  
  for(i in 1:6){
    BetaColumnNames[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ModelLM2", VariablesOfInterest[i], "Beta", sep="_")
  }
  
  print(cbind(colnames(LM2BetasExpression), BetaColumnNames))
  
  colnames(LM2BetasExpression)<-BetaColumnNames
  
  
  LM2DirectionOfEffect<-matrix("Up", nrow(LM2BetasExpression), ncol(LM2BetasExpression))
  
  LM2DirectionOfEffect<-LM2BetasExpression
  LM2DirectionOfEffect[LM2DirectionOfEffect>0]<-"Up"
  LM2DirectionOfEffect[LM2DirectionOfEffect<0]<-"Down"
  #table(LM2DirectionOfEffect[,3])
  #head(cbind(LM2DirectionOfEffect,LM2BetasExpression))
  
  colnames(LM2DirectionOfEffect)
  for(i in 1:6){
    colnames(LM2DirectionOfEffect)[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ModelLM2", VariablesOfInterest[i], "DirectionOfEffect", sep="_")
  }
  
  
  LM2Tstat<-read.table("LM2Tstat.csv", header=T, row.names=1, sep=",")
  colnames(LM2Tstat)
  LM2TstatExpression<-LM2Tstat[,c(2:7)]
  rm(LM2Tstat)
  
  TstatColumnNames<-VariablesOfInterest
  for(i in 1:6){
    TstatColumnNames[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ModelLM2", VariablesOfInterest[i], "Tstat", sep="_")
  }
  
  print(cbind(colnames(LM2TstatExpression), TstatColumnNames))
  
  colnames(LM2TstatExpression)<-TstatColumnNames
  
  LM2pvalues<-read.table("LM2pvaluesRAW.csv", header=T, row.names=1, sep=",")
  
  LM2pvaluesExpression<-LM2pvalues[,c(2:7)]
  rm(LM2pvalues)
  
  PvalueColumnNames<-VariablesOfInterest
  for(i in 1:6){
    PvalueColumnNames[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ModelLM2", VariablesOfInterest[i], "NominalPvalue", sep="_")
  }
  
  print(cbind(colnames(LM2pvaluesExpression), PvalueColumnNames))
  
  colnames(LM2pvaluesExpression)<-PvalueColumnNames
  
  PercentileRank<-LM2pvaluesExpression
  
  for(i in c(1:ncol(LM2pvaluesExpression))){
    PercentileRank[,i]<-rank(as.numeric(LM2pvaluesExpression[,i]), ties.method=("average"))/length(LM2pvaluesExpression[,i])
  }
  
  for(i in 1:6){
    colnames(PercentileRank)[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ModelLM2", VariablesOfInterest[i], "PercentileRank", sep="_")
  }
  
  head(cbind(LM2pvaluesExpression, PercentileRank))
  
  
  #This code may end up crashing - I'm not sure that the files are all named the same for each region. But the calculation of adjusted pvalues is time intensive, so here goes nothing:
  
  #I need to make sure I'm following the same order as VariablesOfInterest
  
  LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
  colnames(LM2testOutputMDDv2)
  head(LM2testOutputMDDv2)
 
  #this is the code that needs to be tweaked depending on what columns of annotation Vikram included:
    ProbeAnnotation<-LM2testOutputMDDv2[, AnnotationColumns]

  colnames(ProbeAnnotation)
  head(ProbeAnnotation)
  
  ProbeAnnotationColNames<-colnames(ProbeAnnotation)
  
  for(i in 1:ncol(ProbeAnnotation)){
    ProbeAnnotationColNames[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ProbeAnnotation", colnames(ProbeAnnotation)[i], sep="_")
  }
  
  print(cbind(ProbeAnnotationColNames, colnames(ProbeAnnotation)))
  
  colnames(ProbeAnnotation)<-ProbeAnnotationColNames
  
  BH_AdjustedPval_FDR<-LM2pvaluesExpression
  
  
  #c("MDD", "Age", "RNAIntegrity", "BrainPH", "Gender", "PMI")
  
  BH_AdjustedPval_FDR[,1]<-LM2testOutputMDDv2[,3]
  rm(LM2testOutputMDDv2)  
  
  
  LM2testOutputAge<-read.table("LM2testOutputAgev2.csv", header=T, row.names=2, sep=",")
  colnames(LM2testOutputAge)
  
  BH_AdjustedPval_FDR[,2]<-LM2testOutputAge[,3]
  rm(LM2testOutputAge)
  
  
  LM2testOutputRNAIntegrity<-read.table("LM2testOutputRNAIntegrityv2.csv", header=T, row.names=2, sep=",")
  colnames(LM2testOutputRNAIntegrity)
  
  BH_AdjustedPval_FDR[,3]<-LM2testOutputRNAIntegrity[,3]
  rm(LM2testOutputRNAIntegrity)
  
  
  LM2testOutputPH<-read.table("LM2testOutputBrainPHv2.csv", header=T, row.names=2, sep=",")
  colnames(LM2testOutputPH)
  
  BH_AdjustedPval_FDR[,4]<-LM2testOutputPH[,3]
  rm(LM2testOutputPH) 
  

  LM2testOutputGender<-read.table("LM2testOutputGenderv2.csv", header=T, row.names=2, sep=",")
  colnames(LM2testOutputGender)
  
  BH_AdjustedPval_FDR[,5]<-LM2testOutputGender[,3]
  rm(LM2testOutputGender)
  
  
  LM2testOutputHoursFinal<-read.table("LM2testOutputHoursFinalv2.csv", header=T, row.names=2, sep=",")
  colnames(LM2testOutputHoursFinal)
  
  BH_AdjustedPval_FDR[,6]<-LM2testOutputHoursFinal[,3]
  rm(LM2testOutputHoursFinal)
  
  colnames(BH_AdjustedPval_FDR)
  
  for(i in 1:6){
    colnames(BH_AdjustedPval_FDR)[i]<-paste(BrainRegion, "Illumina_LCM.AMY_ModelLM2", VariablesOfInterest[i], "FDR", sep="_")
  }
  
  temp<<-data.frame(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
  
}

############End of CODE for function to reformat LCM.AMY results into generalized output######################



#Change to output working directory

setwd("~/Documents/AMY LCM 10nuclei/AAA/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)

ReformatLCM_AMY("AAA", c(12:30, 34:40))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "AAA_Illumina_LCM.AMY_ModelLM2.csv")

rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

#oddly formatted annotation - had to tweak code
setwd("~/Documents/AMY LCM 10nuclei/AB/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("AB", c(12, 14:28, 32:38))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "AB_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)
rm(LM2testOutputMDDv2)

setwd("~/Documents/AMY LCM 10nuclei/AHA/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("AHA", c(12, 14:28, 32:38))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "AHA_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/Basal/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("Basal", c(12, 14:28, 32:39))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "Basal_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/Central/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("Central", c(12:24, 28:34))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "Central_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/CO/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("CO", c(12:31, 34:40))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "CO_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/Lateral/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("Lateral", c(12, 14:28, 31:37))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "Lateral_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/Medial/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("Medial", c(12:26, 28:34))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "Medial_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/PAC/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("PAC", c(12:31, 34:40))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "PAC_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)

setwd("~/Documents/AMY LCM 10nuclei/PL/09 LM2 Output_Best")
LM2testOutputMDDv2<-read.table("LM2testOutputMDDv2.csv", header=T, row.names=2, sep=",")
colnames(LM2testOutputMDDv2)
ReformatLCM_AMY("PL", c(12:31, 34:40))
setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase/Illumina_LCM_AMY")
write.csv(temp, "PL_Illumina_LCM.AMY_ModelLM2.csv")
rm(ProbeAnnotation,LM2BetasExpression, LM2DirectionOfEffect, LM2TstatExpression, LM2pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
rm(temp)
