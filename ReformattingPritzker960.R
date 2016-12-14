#Getting the microarray files ready for the MetaQueryDatabase:

#It would probably be easiest to put all of the probe info, betas, p-values, adj. p-values, t-stats, into a single large matrix. 
#Percentile Rank will need to be calculated, along with a "direction of effect" column, and some adj.p-values.
#File names and column names will need to be standardized.
#Add other forms of gene annotation?


#Let's figure out the code for one Pritzker960 file first, then automate it:

#Pritzker960 data release -> set working directory to brain region-> set working directory to "I am interested in diagnosis"

#Parameters for Analysis (test)
BrainRegion<-"ACG"
DissectionCohorts<-"Yes"

############CODE for function to reformat Pritzker960 results into generalized output######################

ReformatPritzker960<-function(BrainRegion, DissectionCohorts){
  
setwd("./Compiled output")
#list.files()

LM4Betas<-read.table("LM4Betas.csv", header=T, row.names=1, sep=",")
#Hmm... note: I'm not going to be able to automate all of this easily, because only some regions have "EarlyLateCohort" as a variable - although those results may not be worth extracting anyway, right?  Intercept probably isn't worth extracting either...I don't think many of these analyses had pH or Age centered when they were run.

colnames(LM4Betas)
head(LM4Betas)
str(LM4Betas)

LM4BetasExpression<-LM4Betas[,c(2:8)]
colnames(LM4BetasExpression)
rm(LM4Betas)

#I'm a little nervous about automatically renaming these - I should probably print something to the screen so that I double check that it is done properly.

VariablesOfInterest<-c("BP", "MDD", "Schiz", "BrainPH", "Age", "Gender", "PMI")
BetaColumnNames<-VariablesOfInterest

for(i in 1:7){
  BetaColumnNames[i]<-paste(BrainRegion, "Illumina_Pritzker960_ModelLM4", VariablesOfInterest[i], "Beta", sep="_")
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
for(i in 1:7){
  colnames(LM4DirectionOfEffect)[i]<-paste(BrainRegion, "Illumina_Pritzker960_ModelLM4", VariablesOfInterest[i], "DirectionOfEffect", sep="_")
}


LM4Tstat<-read.table("LM4Tstat.csv", header=T, row.names=1, sep=",")
colnames(LM4Tstat)
LM4TstatExpression<-LM4Tstat[,c(2:8)]
rm(LM4Tstat)

TstatColumnNames<-VariablesOfInterest
for(i in 1:7){
 TstatColumnNames[i]<-paste(BrainRegion, "Illumina_Pritzker960_ModelLM4", VariablesOfInterest[i], "Tstat", sep="_")
}

print(cbind(colnames(LM4TstatExpression), TstatColumnNames))

colnames(LM4TstatExpression)<-TstatColumnNames

LM4pvalues<-read.table("LM4pvaluesOutputwithQuality.csv", header=T, row.names=1, sep=",")

#colnames(LM4pvalues)[15]<-"ProbeQuality"

#if(DissectionCohorts=="Yes"){
  #ProbeAnnotation<-LM4pvalues[,c(10:15)] 
#}else{ProbeAnnotation<-LM4pvalues[,c(9:14)] }


LM4pvaluesExpression<-LM4pvalues[,c(2:8)]
rm(LM4pvalues)

PvalueColumnNames<-VariablesOfInterest
for(i in 1:7){
  PvalueColumnNames[i]<-paste(BrainRegion, "Illumina_Pritzker960_ModelLM4", VariablesOfInterest[i], "NominalPvalue", sep="_")
}

print(cbind(colnames(LM4pvaluesExpression), PvalueColumnNames))

colnames(LM4pvaluesExpression)<-PvalueColumnNames

PercentileRank<-LM4pvaluesExpression

for(i in c(1:ncol(LM4pvaluesExpression))){
PercentileRank[,i]<-rank(as.numeric(LM4pvaluesExpression[,i]), ties.method=("average"))/length(LM4pvaluesExpression[,i])
}

for(i in 1:7){
  colnames(PercentileRank)[i]<-paste(BrainRegion, "Illumina_Pritzker960_ModelLM4", VariablesOfInterest[i], "PercentileRank", sep="_")
}

head(cbind(LM4pvaluesExpression, PercentileRank))

getwd()
setwd("../")
getwd()

#This code may end up crashing - I'm not sure that the files are all named the same for each region. But the calculation of adjusted pvalues is time intensive, so here goes nothing:

#I need to make sure I'm following the same order as VariablesOfInterest

LM4testOutputBPv2<-read.table("LM4testOutputBPv2.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputBPv2)

if(DissectionCohorts=="Yes"){
ProbeAnnotation<-LM4testOutputBPv2[, c(10:14, 16:22)]
}else{ProbeAnnotation<-LM4testOutputBPv2[, c(10:13, 15:21)]}

colnames(ProbeAnnotation)
head(ProbeAnnotation)

ProbeAnnotationColNames<-colnames(ProbeAnnotation)

for(i in 1:ncol(ProbeAnnotation)){
  ProbeAnnotationColNames[i]<-paste(BrainRegion, "Illumina_Pritzker960_ProbeAnnotation", colnames(ProbeAnnotation)[i], sep="_")
}

print(cbind(ProbeAnnotationColNames, colnames(ProbeAnnotation)))
  
colnames(ProbeAnnotation)<-ProbeAnnotationColNames

BH_AdjustedPval_FDR<-LM4pvaluesExpression

BH_AdjustedPval_FDR[,1]<-LM4testOutputBPv2[,3]
rm(LM4testOutputBPv2)  

LM4testOutputMDv2<-read.table("LM4testOutputMDv2.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputMDv2)

BH_AdjustedPval_FDR[,2]<-LM4testOutputMDv2[,3]
rm(LM4testOutputMDv2)  

LM4testOutputSchizv2<-read.table("LM4testOutputSchizv2.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputSchizv2)

BH_AdjustedPval_FDR[,3]<-LM4testOutputSchizv2[,3]
rm(LM4testOutputSchizv2) 

LM4testOutputPH<-read.table("LM4testOutputPH.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputPH)

BH_AdjustedPval_FDR[,4]<-LM4testOutputPH[,3]
rm(LM4testOutputPH) 

LM4testOutputAge<-read.table("LM4testOutputAge.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputAge)

BH_AdjustedPval_FDR[,5]<-LM4testOutputAge[,3]
rm(LM4testOutputAge)

LM4testOutputGender<-read.table("LM4testOutputGender.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputGender)

BH_AdjustedPval_FDR[,6]<-LM4testOutputGender[,3]
rm(LM4testOutputGender)

LM4testOutputHoursFinal<-read.table("LM4testOutputHoursFinal.csv", header=T, row.names=2, sep=",")
colnames(LM4testOutputHoursFinal)

BH_AdjustedPval_FDR[,7]<-LM4testOutputHoursFinal[,3]
rm(LM4testOutputHoursFinal)

colnames(BH_AdjustedPval_FDR)

for(i in 1:7){
  colnames(BH_AdjustedPval_FDR)[i]<-paste(BrainRegion, "Illumina_Pritzker960_ModelLM4", VariablesOfInterest[i], "FDR", sep="_")
}

temp<<-data.frame(ProbeAnnotation,LM4BetasExpression, LM4DirectionOfEffect, LM4TstatExpression, LM4pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)

}

############End of CODE for function to reformat Pritzker960 results into generalized output######################


#Change to output working directory
  
write.csv(temp, "ACG_Illumina_Pritzker960_ModelLM4.csv")

rm(ProbeAnnotation,LM4BetasExpression, LM4DirectionOfEffect, LM4TstatExpression, LM4pvaluesExpression, PercentileRank, BH_AdjustedPval_FDR)


ReformatPritzker960("AMY", "Yes")

#Change to output working directory
write.csv(temp, "AMY_Illumina_Pritzker960_ModelLM4.csv")

rm(temp)

ReformatPritzker960("CB", "Yes")

#Change to output working directory
write.csv(temp, "CB_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("DLPFC", "Yes")

#Change to output working directory
write.csv(temp, "DLPFC_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("HC", "Yes")

#Change to output working directory
write.csv(temp, "HC_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("NACC", "Yes")

write.csv(temp, "NACC_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("ANTHAL", "Yes")

write.csv(temp, "ANTHAL_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("MTHAL", "Yes")

write.csv(temp, "MTHAL_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("PCG", "No")
#Crash - I thought there might be some problems there.  
#added some if/else statements to fix it.

write.csv(temp, "PCG_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)

ReformatPritzker960("SCG", "No")

write.csv(temp, "SCG_Illumina_Pritzker960_ModelLM4.csv")
rm(temp)
