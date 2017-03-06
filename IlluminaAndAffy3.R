#creates output for both affy and illumina
#all desired output must be listed
AffyAndIllumina <- function(ILMNDataModel, AFFYDataModel, genesOfInterest, DataSets, AffyDataSets, BrainRegion, variablesofInterest, OutputsInterest){
    
    ILMN <- HiddenDataSet(DataSetModel = ILMNDataModel, geneList = genesOfInterest, datasets = DataSets, brainRegion = BrainRegion, variablesInterest = variablesofInterest, typeOfOutput = OutputsInterest)
    AFFY <- HiddenAffy(AffyDataModel = AFFYDataModel, geneList = genesOfInterest, datasets = DataSets, brainRegion = BrainRegion, variablesInterest = variablesofInterest, typeOfOutput = OutputsInterest)

    JoinedAffyIlmn <- join(AFFY, ILMN, by = "SYMBOLREANNOTATED", type = "full")
    
    output <- list(ILMN, AFFY, JoinedAffyIlmn)
    return (output)
}


#Example Input data

library(plyr)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/thefunction")

illuminaProbeInfo <- read.csv("IlluminaProbeInfo.csv", header = T)
ILMNDataSetModel <- read.csv("DatasetModel.csv", header = T)
AffyDataSetModel <- read.csv("AffyDataSetModel.csv", header = T)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput")
#testgeneList<-as.data.frame(read.csv("TestInput_18_NomP00001.csv", header=T, stringsAsFactors = F))
testgeneList<-as.data.frame(read.csv("TestInput_137_NomP0005.csv", header=T, stringsAsFactors = F))
#testgeneList<-as.data.frame(read.csv("TestInput_714_NomP01.csv", header=T, stringsAsFactors = F))
#testgeneList<-as.data.frame(read.csv("TestInput_1900_nomP05.csv", header=T, stringsAsFactors = F))
#testgeneList<-as.data.frame(read.csv("TestInput_21800.csv", header=T, stringsAsFactors = F))

# setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")
# temp<-read.csv("HC_Affymetrix_Affy6Region_ModelLM4.csv", header=T, stringsAsFactors = F)
# testgeneList<-as.data.frame(temp$SYMBOLREANNOTATED)

testgeneList<-testgeneList[,1]
#colnames(testgeneList)[1]<-"[,1]"
str(testgeneList)


# testgeneList <- t(data.frame( "AAAS", "AACS"))
testgeneList <- c("AAAS", "AACS", "A2M", "A4GNT", "AAMP", "AADAC", "AAK1", "ABAT", "ABCA12", "ABCA2", 
 "ABCA3", "ABCA4", "ABCA5" ,"ABCA6", "ABCA7", "ABCA8", "ABCB1", "ABCB11" ,"ABCB4", "ABCB6", "ABCB7", "ABCB8","EXOSC7", "TES", "TMEM144", "SCUBE2", "UCP2")


testdatasets <- t(data.frame("Pritzker960", "Freeze3", "LCM.AMY", "Affy6Region", "LCM.HC"))
#testBrainRegion <- t(data.frame("HC","DLPFC","NACC","AMY","ACG","CB"))
testBrainRegion <- c("AAA","AB","ACG","AHA","AMY","ANTHAL","Basal","CB","Central","CO","DLPFC","HC","Lateral","Medial","MTHAL","NACC","PAC","PCG","PL","SCG")

testBrainRegion <- c("HC")

testvariablesInterest <- t(data.frame ("MDD", "BP", "SCHIZ", "Age", "BrainPH")) 
#testtypeOfOutput <- t(data.frame("PercentileRank", "Tstat", "NominalPvalue")) #all values must be in one column

testtypeOfOutput <- t(data.frame("Tstat")) #all values must be in one column

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AleksPracticeDatasets/recodeinputgenes")

system.time(
test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = testgeneList, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
)

#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput")

TestOutputAffy<-test[[2]]
TestOutputIllumina<-test[[1]]
TestOutputJoined<-test[[3]]

#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_18_NomP00001")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_137_NomP0005")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_714_NomP01")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_1900_nomP05")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_21800")

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets")
write.csv(TestOutputAffy, "TestOutputAffy.csv")
write.csv(TestOutputIllumina, "TestOutputIllumina.csv")
write.csv(TestOutputJoined, "TestOutputJoined.csv")

#20 brain regions, 21,800 genes, 3 outputs, one variable of interest
# user  system elapsed 
# 141.849   8.216 150.051

testILMN <- HiddenDataSet(DataSetModel = ILMNDataSetModel, geneList = testgeneList, datasets = testdatasets, brainRegion = testBrainRegion, variablesInterest = testvariablesInterest, typeOfOutput = testtypeOfOutput)


testAFFY4 <- HiddenAffy(AffyDataModel = AffyDataSetModel, geneList = testgeneList, datasets = testdatasets, brainRegion = testBrainRegion, variablesInterest = testvariablesInterest, typeOfOutput = testtypeOfOutput)




s <- test[2]
s <- as.data.frame(s)


#############
#Code i found useful
testtype <- c("PercentileRank", "Tstat", "NominalPvalue") #all values must be in one column
unique(grep(paste(testtype,collapse="|"), colnames(testAFFY4), value=F))
#############

#run to add the function to your library
HiddenDataSet <- function(DataSetModel, geneList, datasets, brainRegion, variablesInterest, typeOfOutput){
    # for (j in 1:length(platform)){
    #iterates throught the dataset column of ILMNDataModel
    #pulls out all data related to genes of interest, type of dataset, and brain region
    for (i in 1:length(DataSetModel[,4])){
        if ((DataSetModel[i,4] %in% datasets) & (DataSetModel[i,3] == "Illumina") & (DataSetModel[i,2] %in% brainRegion)){ # taking out for now (DataSetModel[i,4]=="Pritzker960") &
            file <- read.csv(file = toString(DataSetModel[i,1]), stringsAsFactors = FALSE)
            iterationDataSet<-file[as.character(file[,grep("SYMBOLREANNOTATED", colnames(file), perl = T,value = F)]) %in% geneList,]
            
            if(!exists("DataSet")){
                DataSet <- iterationDataSet
            }
            else if(exists("DataSet")){
                DataSet <- join(DataSet, iterationDataSet, type="full", by = "ProbeID")
            }
        }
    }
    #return(DataSet)
    if(!exists("FinalDataSet")){
        FinalDataSet <- DataSet
    }
    #Don't need, will remove and make sure it doesn't affect
    else if(exists("FinalDataSet")){
        FinalDataSet <- join( FinalDataSet, DataSet, type="full" , by = "SYMBOLREANNOTATED")
    }
    #}
    
    #return(FinalDataSet)
    #Grabs all the variable types you want
    for (k in 1:length(variablesInterest)){
        #This first statement adds the Probe id to the dataset so it can be joined with the annotation data
        if (!exists("variableCols")){
            variableCols <- FinalDataSet[,c(grep(variablesInterest[k],c(colnames(FinalDataSet)), perl = T, value = F))]
            variableCols <- cbind(FinalDataSet[,c(grep("ProbeID", c(colnames(FinalDataSet)), perl = T, value = F))], variableCols)
            colnames(variableCols)[1] <- "ProbeID"
        }
        #this only cbinds the next variable to the full Variable dataset
        else if(exists("variableCols")){
            variableColsNext <- FinalDataSet[,c(grep(variablesInterest[k], colnames(FinalDataSet), perl = T, value = F))]
            variableCols <- cbind(variableCols, variableColsNext)
        }
    }
    
    #return (variableCols)
    #Adds In the output types you want
    for (t in 1:length(typeOfOutput)){
        #initializes with probe id to join to the ILMN data
        if (!exists("newerOutput")){
            newerOutput <- variableCols[,c(grep(typeOfOutput[t], c(colnames(variableCols)), perl = T, value = F))]
            newerOutput <- cbind(variableCols[,c(grep("ProbeID", c(colnames(variableCols)), perl = T, value = F))],  newerOutput);
            colnames(newerOutput)[1] <- "ILMN_ProbeID"
            #return (newerOutput)
        }
        #only adds the next ID
        else if (exists("newerOutput")){
            newerOutNext <- variableCols[,c(grep(typeOfOutput[t], c(colnames(variableCols)), perl = T, value = F))]
            newerOutput <- cbind(newerOutput, newerOutNext)
        }
    }
    
    
    
    
    #return (newerOutput)
    
    #joins the annotation
    finalIllumina <- join(illuminaProbeInfo, newerOutput, by = "ILMN_ProbeID", type = "right")
    colnames(finalIllumina)[3] <- "SYMBOLREANNOTATED"
    return (finalIllumina)
}


#run to add AFFY
HiddenAffy <- function(AffyDataModel, geneList, datasets, brainRegion, variablesInterest, typeOfOutput){
    
    for (i in 1:length(AffyDataModel[,4])){
        if ((AffyDataModel[i,4] %in% datasets) & (AffyDataModel[i,3] == "Affymetrix") & (AffyDataModel[i,2] %in% brainRegion)){ # taking out for now (AffyDataModel[i,4]=="Pritzker960") &
            file <- read.csv(file = toString(AffyDataModel[i,1]), stringsAsFactors = FALSE)
            iterationDataSet<-file[as.character(file[,grep("SYMBOLREANNOTATED", colnames(file), perl = T,value = F)]) %in% geneList,]
            
            if(!exists("DataSet")){
                DataSet <- iterationDataSet
            }
            else if(exists("DataSet")){
                DataSet <- join(DataSet, iterationDataSet, type="full", by = "AffyID")
            }
        }
    }
    #return(DataSet)
    if(!exists("FinalDataSet")){
        FinalDataSet <- DataSet
    }
    
    #return (FinalDataSet)
    
    for (k in 1:length(variablesInterest)){
        if (!exists("variableCols")){
            variableCols <- FinalDataSet[,c(grep(variablesInterest[k],c(colnames(FinalDataSet)), perl = T, value = F))]
            variableCols <- cbind(FinalDataSet[,c(grep("AffyID", colnames(FinalDataSet), perl = T, value = F))], variableCols)
            variableCols <- cbind(FinalDataSet[,c(grep("SYMBOLREANNOTATED", colnames(FinalDataSet), perl = T, value = F))], variableCols)
            colnames(variableCols)[grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)] <- "SYMBOLREANNOTATED"
            colnames(variableCols)[grep("AffyID", colnames(variableCols), perl = T, value = F)] <- "AffyID"
            
        }
        else if(exists("variableCols")){
            variableColsNext <- FinalDataSet[,c(grep(variablesInterest[k], colnames(FinalDataSet), perl = T, value = F))]
            variableCols <- cbind(variableCols, variableColsNext)
        }
    }
    #return (variableCols)
    #for (k in 1:length(variablesInterest)){
    #    if (!exists("variableCols")){
    #        variableCols <- FinalDataSet[,c(grep(variablesInterest[k],c(colnames(FinalDataSet)), perl = T, value = F))]
    #        variableCols <- cbind(FinalDataSet[,c(grep("AffyID", colnames(FinalDataSet), perl = T, value = F))], variableCols)
    #        colnames(variableCols)[1] <- "AffyID"
    #    }
    #    else if(exists("variableCols")){
    #        variableColsNext <- FinalDataSet[,c(grep(variablesInterest[k], colnames(FinalDataSet), perl = T, value = F))]
    #        variableCols <- cbind(variableCols, variableColsNext)
    #    }
    #}
    print(str(variableCols))
    #variableCols[,grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)] <- as.character(variableCols[,grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)])
    #print(variableCols)
    
    #for (t in 1:length(typeOfOutput)){
    #    if (!exists("newerOutput")){
    #        newerOutput <- variableCols[,c(grep(toString(typeOfOutput[t]), (colnames(variableCols)), perl = T, value = F))]
    #        newerOutput <- cbind(as.character(variableCols[,c(grep("AffyID", colnames(variableCols), perl = T, value = F))]), newerOutput)
    #        newerOutput <- cbind(as.character(variableCols[,c(grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F))]), newerOutput)
    #        colnames(newerOutput)[3] <- colnames(variableCols)[grep(toString(typeOfOutput[t]), colnames(variableCols), perl = T, value = F)]
    #        
    #        colnames(newerOutput)[1] <- "SYMBOLREANNOTATED"
    #        colnames(newerOutput)[2] <- "AffyID"
    #        #return (newerOutput)
    #    }
    #    else if (exists("newerOutput")){
    #        newerOutNext <- variableCols[,c(grep(toString(typeOfOutput[t]), c(colnames(variableCols)), perl = T, value = F))]
    #        newerOutput <- cbind(newerOutput, newerOutNext)
    #        colnames(newerOutput)[grep("newerOutNext", colnames(newerOutput), perl = T, value = F)] <- colnames(variableCols)[grep(toString(typeOfOutput[t]), colnames(variableCols), perl = T, value = F)]
    #        
    #    }
    #}
    
    newerOutputProbeID <- variableCols[,c(grep("AffyID", colnames(variableCols), perl = T, value = F))]
    newerOutputSYMBOLS <- variableCols[,c(grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F))]
    newerOutPutData <- variableCols[,unique(grep(paste(typeOfOutput,collapse="|"), colnames(variableCols), value=F))]
    newerOutput <- cbind (newerOutputProbeID, newerOutputSYMBOLS, newerOutPutData)
    #return (newerOutput)
    colnames(newerOutput)[grep("newerOutputProbeID", colnames(newerOutput), perl = T, value = F)] <- "AffyID"
    colnames(newerOutput)[grep("newerOutputSYMBOLS", colnames(newerOutput), perl = T, value = F)] <- "SYMBOLREANNOTATED"
    #for (t in 1:length(typeOfOutput)){
    #    if (!exists("newerOutput")){
    #        newerOutput <- variableCols[,c(grep(typeOfOutput[t], c(colnames(variableCols)), perl = T, value = F))]
    #        newerOutput <- cbind(variableCols[,c(grep("AffyID", colnames(variableCols), perl = T, value = F))], newerOutput)
    #        newerOutput <- cbind(variableCols[,c(grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F))], newerOutput)
    #        colnames(newerOutput)[1] <- "SYMBOLREANNOTATED"
    #        colnames(newerOutput)[2] <- "AffyID"
    #    }
    #    else if (exists("newerOutput")){
    #        newerOutNext <- variableCols[,c(grep(typeOfOutput[t], c(colnames(variableCols)), perl = T, value = F))]
    #        newerOutput <- cbind(newerOutput, newerOutNext)
    #    }
    #}
    newerOutput <- as.data.frame(newerOutput)
    return (newerOutput)
}

