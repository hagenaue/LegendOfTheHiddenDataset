
#ReformatLCM_HC

BrainRegion<-"HC"  

  LM5Betas<-read.table("LM5Betas.csv", header=T, row.names=1, sep=",")
  #Hmm... note: I'm not going to be able to automate all of this easily, because only some regions have "EarlyLateCohort" as a variable - although those results may not be worth extracting anyway, right?  Intercept probably isn't worth extracting either...I don't think many of these analyses had pH or Age centered when they were run.
  
  colnames(LM5Betas)
  head(LM5Betas)
  str(LM5Betas)
  
  LM5BetasExpression<-LM5Betas[,c(2:7)]
  colnames(LM5BetasExpression)

  
  #I'm a little nervous about automatically renaming these - I should probably print something to the screen so that I double check that it is done properly.
  
  VariablesOfInterest<-c("BP", "MDD", "Schiz","BrainPH", "Age", "PMI")
  BetaColumnNames<-VariablesOfInterest
  
  for(i in 1:6){
    BetaColumnNames[i]<-paste(BrainRegion, "Illumina_LCM.HC_ModelLM5", VariablesOfInterest[i], "Beta", sep="_")
  }
  
  print(cbind(colnames(LM5BetasExpression), BetaColumnNames))
  
  colnames(LM5BetasExpression)<-BetaColumnNames
  
  
  LM5DirectionOfEffect<-matrix("Up", nrow(LM5BetasExpression), ncol(LM5BetasExpression))
  
  LM5DirectionOfEffect<-LM5BetasExpression
  LM5DirectionOfEffect[LM5DirectionOfEffect>0]<-"Up"
  LM5DirectionOfEffect[LM5DirectionOfEffect<0]<-"Down"
  #table(LM5DirectionOfEffect[,3])
  #head(cbind(LM5DirectionOfEffect,LM5BetasExpression))
  
  colnames(LM5DirectionOfEffect)
  for(i in 1:6){
    colnames(LM5DirectionOfEffect)[i]<-paste(BrainRegion, "Illumina_LCM.HC_ModelLM5", VariablesOfInterest[i], "DirectionOfEffect", sep="_")
  }
  
  
  #LM5Tstat<-read.table("LM5Tstat.csv", header=T, row.names=1, sep=",")
  colnames(LM5Tstat)
  LM5TstatExpression<-LM5Tstat[,c(2:7)]
  
  TstatColumnNames<-VariablesOfInterest
  for(i in 1:6){
    TstatColumnNames[i]<-paste(BrainRegion, "Illumina_LCM.HC_ModelLM5", VariablesOfInterest[i], "Tstat", sep="_")
  }
  
  print(cbind(colnames(LM5TstatExpression), TstatColumnNames))
  
  colnames(LM5TstatExpression)<-TstatColumnNames
  
  #LM5pvalues<-read.table("LM5pvaluesRAW.csv", header=T, row.names=1, sep=",")
  
  LM5pvaluesExpression<-LM5pvalues[,c(2:7)]
  #rm(LM5pvalues)
  
  PvalueColumnNames<-VariablesOfInterest
  for(i in 1:6){
    PvalueColumnNames[i]<-paste(BrainRegion, "Illumina_LCM.HC_ModelLM5", VariablesOfInterest[i], "NominalPvalue", sep="_")
  }
  
  print(cbind(colnames(LM5pvaluesExpression), PvalueColumnNames))
  
  colnames(LM5pvaluesExpression)<-PvalueColumnNames
  
  PercentileRank<-LM5pvaluesExpression
  
  for(i in c(1:ncol(LM5pvaluesExpression))){
    PercentileRank[,i]<-rank(as.numeric(LM5pvaluesExpression[,i]), ties.method=("average"))/length(LM5pvaluesExpression[,i])
  }
  
  for(i in 1:6){
    colnames(PercentileRank)[i]<-paste(BrainRegion, "Illumina_LCM.HC_ModelLM5", VariablesOfInterest[i], "PercentileRank", sep="_")
  }
  
  head(cbind(LM5pvaluesExpression, PercentileRank))
  
  
  #This code may end up crashing - I'm not sure that the files are all named the same for each region. But the calculation of adjusted pvalues is time intensive, so here goes nothing:
  
  #I need to make sure I'm following the same order as VariablesOfInterest
  #[1] "BP"      "MDD"     "Schiz"   "BrainPH" "Age"     "PMI"
  
  
  LM5testOutputBPv2<-read.table("LM5testOutputBPv2   THIS ONE.csv", header=T, row.names=2, sep=",")
  colnames(LM5testOutputBPv2)
  
    ProbeAnnotation<-LM5testOutputBPv2[, c(12:24)]

  colnames(ProbeAnnotation)
  head(ProbeAnnotation)
  
  ProbeAnnotationColNames<-colnames(ProbeAnnotation)
  
  for(i in 1:ncol(ProbeAnnotation)){
    ProbeAnnotationColNames[i]<-paste(BrainRegion, "Illumina_LCM.HC_ProbeAnnotation", colnames(ProbeAnnotation)[i], sep="_")
  }
  
  print(cbind(ProbeAnnotationColNames, colnames(ProbeAnnotation)))
  
  colnames(ProbeAnnotation)<-ProbeAnnotationColNames
  
  BH_AdjustedPval_FDR<-LM5pvaluesExpression
  
  BH_AdjustedPval_FDR[,1]<-LM5testOutputBPv2[,3]
  rm(LM5testOutputBPv2)  
  
  
  LM5testOutputMDDv2<-read.table("LM5testOutputMDDv2 THIS ONE.csv", header=T, row.names=2, sep=",")
  colnames(LM5testOutputMDDv2)
  head(LM5testOutputMDDv2)

  BH_AdjustedPval_FDR[,2]<-LM5testOutputMDDv2[,3]
  rm(LM5testOutputMDDv2)  
  
  
  LM5testOutputSchizv2<-read.table("LM5testOutputSchizv2   THIS ONE.csv", header=T, row.names=2, sep=",")
  colnames(LM5testOutputSchizv2)
  head(LM5testOutputSchizv2)
  
  BH_AdjustedPval_FDR[,3]<-LM5testOutputSchizv2[,3]
  rm(LM5testOutputSchizv2)  
  
  
  LM5testOutputPH<-read.table("LM5testOutputPHv2.csv", header=T, row.names=2, sep=",")
  colnames(LM5testOutputPH)
  
  BH_AdjustedPval_FDR[,4]<-LM5testOutputPH[,3]
  rm(LM5testOutputPH) 
  
  
  LM5testOutputAge<-read.table("LM5testOutputAgev2.csv", header=T, row.names=2, sep=",")
  colnames(LM5testOutputAge)
  
  BH_AdjustedPval_FDR[,5]<-LM5testOutputAge[,3]
  rm(LM5testOutputAge)
  
  
  LM5testOutputHoursFinal<-read.table("LM5testOutputHoursFinalv2.csv", header=T, row.names=2, sep=",")
  colnames(LM5testOutputHoursFinal)
  
  BH_AdjustedPval_FDR[,6]<-LM5testOutputHoursFinal[,3]
  rm(LM5testOutputHoursFinal)
  
  colnames(BH_AdjustedPval_FDR)
  
  
  for(i in 1:6){
    colnames(BH_AdjustedPval_FDR)[i]<-paste(BrainRegion, "Illumina_LCM.HC_ModelLM5", VariablesOfInterest[i], "FDR", sep="_")
  }
  
  temp<<-data.frame(ProbeAnnotation,LM5BetasExpression, LM5DirectionOfEffect, LM5TstatExpression, LM5pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)
  


############End of CODE for function to reformat LCM.HC results into generalized output######################
  write.csv(temp, "HC_Illumina_LCM.HC_ModelLM5.csv")
  
  setwd("/Users/mhh/Documents/R Code/MakingAMetaQueryDatabase")
  write.csv(temp, "HC_Illumina_LCM.AMY_ModelLM5.csv")
  
  
  