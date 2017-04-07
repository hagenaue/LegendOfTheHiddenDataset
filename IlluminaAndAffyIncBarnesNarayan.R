#creates output for both affy and illumina
#Parent function, creates a list of three datasets in the order: illumina, affy, illumina+Affy by genesymbol.
AffyAndIllumina <- function(ILMNDataModel, AFFYDataModel, genesOfInterest, DataSets, BrainRegion, variablesofInterest, OutputsInterest){
    
    if("LCM.AMY" %in% DataSets | "Freeze3" %in% DataSets | "Pritzker960" %in% DataSets | "LCM.HC"%in% DataSets)
        ILMN <- HiddenIllumina(DataSetModel = ILMNDataModel, geneList = genesOfInterest, datasets = DataSets, brainRegion = BrainRegion, variablesInterest = variablesofInterest, typeOfOutput = OutputsInterest)
    if("Affy6Region" %in% DataSets)
        AFFY <- HiddenAffy(AffyDataModel = AFFYDataModel, geneList = genesOfInterest, datasets = DataSets, brainRegion = BrainRegion, variablesInterest = variablesofInterest, typeOfOutput = OutputsInterest)
    if("Barnes" %in% DataSets | "Narayan" %in% DataSets)
        narayanbarnes <- BarnesNarayan(AffyDataModel = AFFYDataModel, geneList = genesOfInterest, datasets = DataSets, brainRegion = BrainRegion, variablesInterest = variablesofInterest,typeOfOutput = OutputsInterest)
    
    if (exists("AFFY") & exists("narayanbarnes") & exists("ILMN")){
        joinedAffy <- join(AFFY, narayanbarnes, by = "SYMBOLREANNOTATED", type = "full")
        JoinedAffyIlmn <- join(joinedAffy, ILMN, by = "SYMBOLREANNOTATED", type = "full")
        output <- list(ILMN, joinedAffy, JoinedAffyIlmn)
        return (output)
        break
    }
    if (exists("AFFY") & exists("ILMN")){
        joinedAffyIlmn <- join(AFFY, ILMN, by = "SYMBOLREANNOTATED", type = "full")
        output <- list(ILMN, AFFY, joinedAffyIlmn)
        return (output)
        break
    }
    if (exists("joinedAffy") & !exists("ILMN")){
        return (joinedAffy)
        break
    }
    if(exists("AFFY") & !exists("joinedAffy")){
        return (AFFY)
        break
    }
    if(exists("ILMN") & !exists("AFFY"))
        return (ILMN)
}
#run to add the function to your library
HiddenIllumina <- function(DataSetModel, geneList, datasets, brainRegion, variablesInterest, typeOfOutput){
    #iterates throught the dataset column of ILMNDataModel
    #pulls out all data related to genes of interest, type of dataset, and brain region
    for (i in 1:length(DataSetModel[,4])){
        if ((DataSetModel[i,4] %in% datasets) & (DataSetModel[i,3] == "Illumina") & (DataSetModel[i,2] %in% brainRegion)){ # taking out for now (DataSetModel[i,4]=="Pritzker960") &
            file <- read.csv(file = toString(DataSetModel[i,1]), stringsAsFactors = FALSE)
            iterationDataSet<-file[as.character(file[,grep("SYMBOLREANNOTATED", colnames(file), ignore.case = TRUE, perl = T,value = F)]) %in% geneList,]
            
            if(!exists("DataSet")){
                DataSet <- iterationDataSet
            }
            else if(exists("DataSet")){
                DataSet <- join(DataSet, iterationDataSet, type="full", by = "ProbeID")
            }
        }
    }
    FinalDataSet <- DataSet
    
    #Grabs all the variable types you want
    for (k in 1:length(variablesInterest)){
        #This first statement adds the Probe id to the dataset so it can be joined with the annotation data
        if (!exists("variableCols")){
            variableCols <- FinalDataSet[,c(grep(variablesInterest[k],c(colnames(FinalDataSet)), ignore.case = TRUE, perl = T, value = F))]
            variableCols <- cbind(FinalDataSet[,c(grep("ProbeID", c(colnames(FinalDataSet)), ignore.case = TRUE, perl = T, value = F))], variableCols)
            colnames(variableCols)[1] <- "ProbeID"
        }
        #this only cbinds the next variable columns to the full Variable dataset
        else if(exists("variableCols")){
            variableColsNext <- FinalDataSet[,c(grep(variablesInterest[k], colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))]
            variableCols <- cbind(variableCols, variableColsNext)
        }
    }
    #Adds In the output types you want
    for (t in 1:length(typeOfOutput)){
        #initializes with probe id to join to the ILMN data and the first output type
        if (!exists("newerOutput")){
            newerOutput <- variableCols[,c(grep(typeOfOutput[t], c(colnames(variableCols)), ignore.case = TRUE, perl = T, value = F))]
            newerOutput <- cbind(variableCols[,c(grep("ProbeID", c(colnames(variableCols)), ignore.case = TRUE, perl = T, value = F))],  newerOutput);
            colnames(newerOutput)[1] <- "ILMN_ProbeID"
            #return (newerOutput)
        }
        #only adds the next type of output
        else if (exists("newerOutput")){
            newerOutNext <- variableCols[,c(grep(typeOfOutput[t], c(colnames(variableCols)), ignore.case = TRUE, perl = T, value = F))]
            newerOutput <- cbind(newerOutput, newerOutNext)
        }
    }
    
    #joins the annotation and data together
    finalIllumina <- join(illuminaProbeInfo, newerOutput, by = "ILMN_ProbeID", type = "right")
    colnames(finalIllumina)[3] <- "SYMBOLREANNOTATED"
    return (finalIllumina)
}
#run to add AFFY
HiddenAffy <- function(AffyDataModel, geneList, datasets, brainRegion, variablesInterest, typeOfOutput){
    #Grabs all data related to dataset input, brainRegion, and input genes, from each file
    for (i in 1:length(AffyDataModel[,4])){
        if ((AffyDataModel[i,4] == "Affy6Region") & (AffyDataModel[i,3] == "Affymetrix") & (AffyDataModel[i,2] %in% brainRegion)){ # taking out for now (AffyDataModel[i,4]=="Pritzker960") &
            #initializes the dataset with genesymbol column and the first file 
            file <- read.csv(file = toString(AffyDataModel[i,1]), stringsAsFactors = FALSE)
            iterationDataSet<-file[as.character(file[,grep("SYMBOLREANNOTATED", colnames(file), ignore.case = TRUE, perl = T,value = F)]) %in% geneList,]
            
            #if already initialized, just reads next file rather than adding symbol column
            if(!exists("DataSet")){
                DataSet <- iterationDataSet
            }
            #joins the two datasets by id
            else if(exists("DataSet")){
                DataSet <- join(DataSet, iterationDataSet, type="full", by = "AffyID")
            }
        }
    }
    
    FinalDataSet <- DataSet
    
    #pulls out each individual variable of interest and binds them together
    for (k in 1:length(variablesInterest)){
        #initializes the new data with ID and Symbol columns
        if (!exists("variableCols")){
            variableCols <- FinalDataSet[,c(grep(variablesInterest[k],c(colnames(FinalDataSet)), ignore.case = TRUE, perl = T, value = F))]
            variableCols <- cbind(FinalDataSet[,c(grep("AffyID", colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))], variableCols)
            variableCols <- cbind(FinalDataSet[,c(grep("SYMBOLREANNOTATED", colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))], variableCols)
            colnames(variableCols)[grep("SYMBOLREANNOTATED", colnames(variableCols), ignore.case = TRUE, perl = T, value = F)] <- "SYMBOLREANNOTATED"
            colnames(variableCols)[grep("AffyID", colnames(variableCols), ignore.case = TRUE, perl = T, value = F)] <- "AffyID"
            
        }
        #adds the next variable to the end of the dataset
        else if(exists("variableCols")){
            variableColsNext <- FinalDataSet[,c(grep(variablesInterest[k], colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))]
            variableCols <- cbind(variableCols, variableColsNext)
        }
    }
    #return (variableCols)
    #Test Stuff #######
    #return (variableCols)
    #print(str(variableCols))
    #variableCols[,grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)] <- as.character(variableCols[,grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)])
    #print(variableCols)
    #return(variableCols)
    #separates the probeID, Gene symbol, and data columns.
    newerOutputProbeID <- variableCols[,c(grep("AffyID", colnames(variableCols), ignore.case = TRUE, perl = T, value = F))]
    newerOutputSYMBOLS <- variableCols[,c(grep("SYMBOLREANNOTATED", colnames(variableCols), ignore.case = TRUE, perl = T, value = F))]
    #below pulls the output data you want and puts it into the final dataset with gene symbol and probe column
    newerOutPutData <- variableCols[,unique(grep(paste(typeOfOutput,collapse="|"), colnames(variableCols), ignore.case = TRUE, value=F))]
    #return (newerOutPutData)
    #Binds the probe id, gene symbol and data columns together

    newerOutput <-  data.frame(newerOutputProbeID, newerOutputSYMBOLS, newerOutPutData)
    

    colnames(newerOutput)[grep("newerOutputProbeID", colnames(newerOutput), ignore.case = TRUE, perl = T, value = F)] <- "AffyID"
    colnames(newerOutput)[grep("newerOutputSYMBOLS", colnames(newerOutput), ignore.case = TRUE, perl = T, value = F)] <- "SYMBOLREANNOTATED"
    
    if(is.null(colnames(newerOutPutData))){
        colnames(newerOutput)[3] <- colnames(variableCols[unique(grep(paste(typeOfOutput,collapse="|"), colnames(variableCols), ignore.case = TRUE, value=F))])
    }
    newerOutput <- as.data.frame(newerOutput)
    return (newerOutput)
}
BarnesNarayan <- function(AffyDataModel, geneList, datasets, brainRegion, variablesInterest, typeOfOutput){   
    #Grabs all data related to dataset input, brainRegion, and input genes, from each file
    for (i in 1:length(AffyDataModel[,4])){
        if ((AffyDataModel[i,4] %in% datasets) & (AffyDataModel[i,3] == "Affymetrix") & (AffyDataModel[i,2] %in% brainRegion) & (AffyDataModel[i,4] == "Barnes" | AffyDataModel[i,4] =="Narayan")){ 
            #initializes the dataset with genesymbol column and the first file 
            file <- read.csv(file = toString(AffyDataModel[i,1]), stringsAsFactors = FALSE)
            iterationDataSet<-file[as.character(file[,grep("SYMBOLREANNOTATED", colnames(file), ignore.case = TRUE, perl = T,value = F)]) %in% geneList,]
            
            #if already initialized, just reads next file rather than adding symbol column
            if(!exists("DataSet")){
                DataSet <- iterationDataSet
            }
            #joins the two datasets by id
            else if(exists("DataSet")){
                DataSet <- join(DataSet, iterationDataSet, type="full", by = "SYMBOLREANNOTATED")
            }
        }
    }

    FinalDataSet <- DataSet
    #return (FinalDataSet)
    #pulls out each individual variable of interest and binds them together
    for (k in 1:length(variablesInterest)){
        #initializes the new data with ID and Symbol columns
        if (!exists("variableCols")){
            variableCols <- FinalDataSet[,c(grep(variablesInterest[k],c(colnames(FinalDataSet)), ignore.case = TRUE, perl = T, value = F))]
            variableCols <- cbind(FinalDataSet[,c(grep("EntrezGeneID", colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))], variableCols)
            variableCols <- cbind(FinalDataSet[,c(grep("ProbesetID", colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))], variableCols)
            variableCols <- cbind(FinalDataSet[,c(grep("SYMBOLREANNOTATED", colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))], variableCols)
            colnames(variableCols)[grep("SYMBOLREANNOTATED", colnames(variableCols), ignore.case = TRUE, perl = T, value = F)] <- "SYMBOLREANNOTATED"
            colnames(variableCols)[grep("EntrezGeneID", colnames(variableCols), ignore.case = TRUE, perl = T, value = F)] <- "EntrezGeneID"
            colnames(variableCols)[grep("ProbesetID", colnames(variableCols), ignore.case = TRUE, perl = T, value = F)] <- "ProbesetID"
        }
        #adds the next variable to the end of the dataset
        else if(exists("variableCols")){
            variableColsNext <- FinalDataSet[,c(grep(variablesInterest[k], colnames(FinalDataSet), ignore.case = TRUE, perl = T, value = F))]
            variableCols <- cbind(variableCols, variableColsNext)
        }
    }
    
    #Test Stuff #######
    #return (variableCols)
    #print(str(variableCols))
    #variableCols[,grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)] <- as.character(variableCols[,grep("SYMBOLREANNOTATED", colnames(variableCols), perl = T, value = F)])
    #print(variableCols)
    
    #separates the probeID, Gene symbol, and data columns.
    newerOutputProbesetID <- variableCols[,c(grep("ProbesetID", colnames(variableCols), ignore.case = TRUE, perl = T, value = F))]
    newerOutputEntrezID <- variableCols[,c(grep("EntrezGeneID", colnames(variableCols), ignore.case = TRUE, perl = T, value = F))]
    newerOutputSYMBOLS <- variableCols[,c(grep("SYMBOLREANNOTATED", colnames(variableCols), ignore.case = TRUE, perl = T, value = F))]
    #below pulls the output data you want and puts it into the final dataset with gene symbol and probe column
    newerOutPutData <- variableCols[,unique(grep(paste(typeOfOutput,collapse="|"), colnames(variableCols), ignore.case = TRUE, value=F))]
    
    #Binds the probe id, gene symbol and data columns together
    newerOutput <-  data.frame(newerOutputEntrezID, newerOutputProbesetID, newerOutputSYMBOLS, newerOutPutData)
    colnames(newerOutput)[grep("newerOutputProbesetID", colnames(newerOutput), ignore.case = TRUE, perl = T, value = F)] <- "ProbesetID"
    colnames(newerOutput)[grep("newerOutputEntrezID", colnames(newerOutput), ignore.case = TRUE, perl = T, value = F)] <- "EntrezID"
    colnames(newerOutput)[grep("newerOutputSYMBOLS", colnames(newerOutput), ignore.case = TRUE, perl = T, value = F)] <- "SYMBOLREANNOTATED"
    
    newerOutput <- as.data.frame(newerOutput)
    return (newerOutput)
}
#############################################3
#Example Function usage:
library(plyr)
#3 necessary info files
illuminaProbeInfo <- read.csv("IlluminaProbeInfo.csv", header = T)
ILMNDataSetModel <- read.csv("DatasetModel.csv", header = T)
AffyDataSetModel <- read.csv("AffyDataSetModel.csv", header = T)
#gene list, datasets, BrainRegions, Variables, Types Output Declared
#can be one column of a dataframe or a cbound list
testgeneList <- c("AAAS", "AACS", "A2M", "A4GNT", "AAMP", "AADAC", "AAK1", "ABAT", "ABCA12", "ABCA2", "ABCA3", "ABCA4")
testdatasets <- c("Pritzker960", "Freeze3", "LCM.AMY", "Affy6Region", "Barnes", "Narayan")
testBrainRegion <- c("AAA","AB","ACG","AHA","AMY","DLPFC", "TemporalCortex")
testvariablesInterest <- c("MDD")
testtypeOfOutput <- c("PercentileRank", "Tstat", "NominalPvalue")
#Function Call
test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = testgeneList, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)

ILMN <- test[2]
ILMN <- as.data.frame(ILMN)
#################################################################
#ABSOLUTELY REQUIRED DATA BELOW RUN BEFORE TRYING TO USE FUNCTION
#make sure all data is within your working directory
library(plyr)

illuminaProbeInfo <- read.csv("IlluminaProbeInfo.csv", header = T)
ILMNDataSetModel <- read.csv("DatasetModel.csv", header = T)
AffyDataSetModel <- read.csv("AffyDataSetModel.csv", header = T)

#################################################################
#Input gene list must be a vector of characters. Tnis can come from a file column but you must convert it to a vector
setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput")
###############reads in dataframe of genes######################
#testgeneList<-as.data.frame(read.csv("TestInput_18_NomP00001.csv", header=T, stringsAsFactors = F))
#testgeneList<-as.data.frame(read.csv("TestInput_137_NomP0005.csv", header=T, stringsAsFactors = F))
#testgeneList<-as.data.frame(read.csv("TestInput_714_NomP01.csv", header=T, stringsAsFactors = F))
#testgeneList<-as.data.frame(read.csv("TestInput_1900_nomP05.csv", header=T, stringsAsFactors = F))
testgeneList<-as.data.frame(read.csv("TestInput_21800.csv", header=T, stringsAsFactors = F))
##############Converts DataFrame to vector######################
testgeneList<-testgeneList[,1]
str(testgeneList)

#Can also create a concatonated list of chars
#testgeneList <- c("AAAS", "AACS", "A2M", "A4GNT", "AAMP", "AADAC", "AAK1", "ABAT", "ABCA12", "ABCA2", "ABCA3", "ABCA4", "ABCA5" ,"ABCA6", "ABCA7", "ABCA8", "ABCB1", "ABCB11" ,"ABCB4", "ABCB6", "ABCB7", "ABCB8","EXOSC7", "TES", "TMEM144", "SCUBE2", "UCP2")
#################################################################

#Dataset list Can be one column of a dataframe or a concatonated list
testdatasets <- t(data.frame("Pritzker960", "Freeze3", "LCM.AMY", "Affy6Region"))
testdatasets <- c("Pritzker960", "Freeze3", "LCM.AMY", "Affy6Region")

#These are all 20 brain regions. Any region you want must be listed, anything not listed will be left out
testBrainRegion <- c("AAA","AB","ACG","AHA","AMY","ANTHAL","Basal","CB","Central","CO","DLPFC","HC","Lateral","Medial","MTHAL","NACC","PAC","PCG","PL","SCG")


#variables can be either a list or column of dataframe
testvariablesInterest <- t(data.frame ("MDD")) 
testvariablesInterest <- c("MDD")
#List of all possible variables : MDD, Age, RNAIntegrity, BrainPH, Gender, PMI, BP, SCHIZ, AgonalFactor

testtypeOfOutput <- t(data.frame("PercentileRank", "Tstat", "NominalPvalue")) #all values must be in one column
testtypeOfOutput <- c("PercentileRank", "Tstat", "NominalPvalue")
#List of all possible outputs: Beta, DirectionOfEffect, Tstat, NominalPvalue, PercentileRank, FDR


setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AleksPracticeDatasets/recodeinputgenes")

system.time(
    test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = testgeneList, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
)


#OUTPUTTING THE FILE
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput")

#Accessing individual datasets
TestOutputAffy<-test[[2]]
TestOutputIllumina<-test[[1]]
TestOutputJoined<-test[[3]]

#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_18_NomP00001")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_137_NomP0005")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_714_NomP01")
#setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_1900_nomP05")
setwd("~/Documents/R Code/MakingAMetaQueryDatabase/TestOutput/TestInput_21800")


write.csv(TestOutputAffy, "TestOutputAffy.csv")
write.csv(TestOutputIllumina, "TestOutputIllumina.csv")
write.csv(TestOutputJoined, "TestOutputJoined.csv")

#20 brain regions, 21,800 genes, 3 outputs, one variable of interest
# user  system elapsed 
# 141.849   8.216 150.051


#Calls to hiddenDataSet and hiddenAffy separate from the joined file
testILMN <- HiddenIllumina(DataSetModel = ILMNDataSetModel, geneList = testgeneList, datasets = testdatasets, brainRegion = testBrainRegion, variablesInterest = testvariablesInterest, typeOfOutput = testtypeOfOutput)
testAFFY4 <- HiddenAffy(AffyDataModel = AffyDataSetModel, geneList = testgeneList, datasets = testdatasets, brainRegion = testBrainRegion, variablesInterest = testvariablesInterest, typeOfOutput = testtypeOfOutput)





#############
#Code I found useful
testtype <- c("PercentileRank", "Tstat", "NominalPvalue") #all values must be in one column
unique(grep(paste(testtype,collapse="|"), colnames(testAFFY4), value=F))
#############
