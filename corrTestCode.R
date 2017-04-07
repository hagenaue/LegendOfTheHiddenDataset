illuminaProbeInfo <- read.csv("IlluminaProbeInfo.csv", header = T)
ILMNDataSetModel <- read.csv("DatasetModel.csv", header = T)
AffyDataSetModel <- read.csv("AffyDataSetModel.csv", header = T)

testdatasets <- c("Pritzker960", "Freeze3", "LCM.AMY", "LCM.HC", "Affy6Region")
testtypeOfOutput <- c("Tstat")
allVariables <- c("MDD", "Age", "Gender", "BP", "SCHIZ")
allRegions <- c("DLPFC", "HC", "NACC", "ACG", "CB", "AMY")
dir.create("corrFolder")
for (i in 1:length(allRegions)){
    for (j in 1:length(allVariables)){
        test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = toString(allRegions[i]), variablesofInterest = toString(allVariables[j]), OutputsInterest = testtypeOfOutput)
        data <- test[3]
        data <- as.data.frame(data)
        tStatCols <- data[,grep("Tstat", colnames(data), perl = T, ignore.case = T, value = F)]
        tStatCorMatrix <- cor(tStatCols, use = "pairwise.complete.obs")
        title <- c(toString(allRegions[i]), toString(allVariables[j]), "corr.csv")
        setwd("corrFolder")
        write.csv(tStatCorMatrix, toString(title))
        setwd("..")
    }
}

library(plyr)

#Hippocampal age correlation
testBrainRegion <- c("HC")
testvariablesInterest <- c("Age")
test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
data <- test[3]
data <- as.data.frame(data)
tStatCols <- data[,grep("Tstat", colnames(data), perl = T, ignore.case = T, value = F)]
tStatCorMatrix <- cor(tStatCols, use = "complete.obs")

#Hippocampal MDD corr
testBrainRegion <- c("HC")
testvariablesInterest <- c("MDD")
test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = "AAA", variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
data <- test[3]
data <- as.data.frame(data)
tStatCols <- data[,grep("Tstat", colnames(data), perl = T, ignore.case = T, value = F)]
tStatCorMatrix <- cor(tStatCols, use = "complete.obs")

#Hippocampal BP corr
testBrainRegion <- c("HC")
testvariablesInterest <- c("BP")
test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
data <- test[3]
data <- as.data.frame(data)
tStatCols <- data[,grep("Tstat", colnames(data), perl = T, ignore.case = T, value = F)]
tStatCorMatrix <- cor(tStatCols, use = "complete.obs")





dataILMN <- as.data.frame(test[1])

dataAffy <- as.data.frame(test[2])

testAFFY3 <- HiddenAffy(AffyDataModel = AffyDataSetModel, geneList = genes, datasets = testdatasets, brainRegion = testBrainRegion, variablesInterest = testvariablesInterest, typeOfOutput = testtypeOfOutput)

test <- as.data.frame(testAFFY2)



































genes <- c("A1CF",
"A2M",
"A4GALT",
"A4GNT",
"AAAS",
"AACS",
"AADAC",
"AAGAB",
"AAK1",
"AAMP",
"AANAT",
"AARS",
"AARSD1",
"AASDHPPT",
"AASS",
"AATF",
"AATK",
"ABAT",
"ABCA1",
"ABCA12",
"ABCA2",
"ABCA3",
"ABCA4",
"ABCA5",
"ABCA6",
"ABCA7",
"ABCA8",
"ABCB1",
"ABCB11",
"ABCB4",
"ABCB6",
"ABCB7",
"ABCB8",
"ABCB9",
"ABCC1",
"ABCC10",
"ABCC2",
"ABCC3",
"ABCC4",
"ABCC5",
"ABCC6",
"ABCC8",
"ABCC9",
"ABCD1",
"ABCD2",
"ABCD3",
"ABCD4",
"ABCE1",
"ABCF1",
"ABCF2",
"ABCF3",
"ABCG1",
"ABCG2",
"ABCG4",
"ABCG5",
"ABHD10",
"ABHD11",
"ABHD14A",
"ABHD2",
"ABHD3",
"ABHD4",
"ABHD5",
"ABHD6",
"ABHD8",
"ABI1",
"ABI2",
"ABL1",
"ABL2",
"ABLIM1",
"ABLIM3",
"ABO")
