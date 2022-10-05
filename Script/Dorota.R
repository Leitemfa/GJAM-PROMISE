#Importing data
library(readr)
Phenotypes_2wpi_Cluefield <- read_csv("Dorota/20210326-Phenotypes-2wpi-Cluefield.csv")
View(Phenotypes_2wpi_Cluefield)

Phenotypes_3wpi_Cluefield <- read_csv("Dorota/20210326-Phenotypes-3wpi-Cluefield.csv")
View(Phenotypes_2wpi_Cluefield)

#Bacteria
OTU_Table_Bacteria_Cluefield_bulksoil_subset <- read_csv("OTU_Table_Bacteria_Cluefield_bulksoil_subset.csv")
# View(OTU_Table_Bacteria_Cluefield_bulksoil_subset)

OTU_Table_Bacteria_Cluefield_Endophytes_root_soil_subset <- read_csv("OTU_Table_Bacteria_Cluefield_Endophytes_root_soil_subset.csv")
# View(OTU_Table_Bacteria_Cluefield_Endophytes_root_soil_subset)

OTU_Table_Bacteria_Cluefield_Endophytes_root_sand_subset <- read_csv("OTU_Table_Bacteria_Cluefield_Endophytes_root_sand_subset.csv")
# View(OTU_Table_Bacteria_Cluefield_Endophytes_root_soil_subset)

OTU_Table_Bacteria_Cluefield_Rhizosphere_subset <- read_csv("OTU_Table_Bacteria_Cluefield_Rhizosphere_subset.csv")
# View(OTU_Table_Bacteria_Cluefield_Rhizosphere_subset)

#Fungi
OTU_Table_Fungi_Cluefield_bulksoil_subset <- read_csv("OTU_Table_Fungi_Cluefield_bulksoil_subset.csv")
# View(OTU_Table_Fungi_Cluefield_bulksoil_subset)

OTU_Table_Fungi_Cluefield_Endophytes_root_soil_subset <- read_csv("OTU_Table_Fungi_Cluefield_Endophytes_root_soil-subset.csv")
# View(OTU_Table_Fungi_Cluefield_Endophytes_root_soil_subset)

OTU_Table_Fungi_Cluefield_Endophytes_root_sand_subset <- read_csv("OTU_Table_Fungi_Cluefield_Endophytes_root_sand-subset.csv")
# View(OTU_Table_Fungi_Cluefield_Endophytes_root_soil_subset)

OTU_Table_Fungi_Cluefield_Rhizosphere_subset <- read_csv("OTU_Table_Fungi_Cluefield_Rhizospere-subset.csv")
# View(OTU_Table_Fungi_Cluefield_Rhizosphere_subset)

#Data organization
library(tidyr)
library(dplyr)
TaxInfo <- OTU_Table_Bacteria_Cluefield_bulksoil_subset[,c(1,50)]
TaxInfoF <- OTU_Table_Fungi_Cluefield_bulksoil_subset[,c(1,50)]
#input data
otu_data <- OTU_Table_Bacteria_Cluefield_bulksoil_subset
otu_data <- OTU_Table_Bacteria_Cluefield_Endophytes_root_soil_subset
otu_data <- OTU_Table_Bacteria_Cluefield_Endophytes_root_sand_subset %>% 
  left_join(TaxInfo, by = "#OTU ID")
otu_data <- OTU_Table_Bacteria_Cluefield_Rhizosphere_subset

otu_data <- OTU_Table_Fungi_Cluefield_bulksoil_subset
otu_data <- OTU_Table_Fungi_Cluefield_Endophytes_root_soil_subset
otu_data <- OTU_Table_Fungi_Cluefield_Endophytes_root_sand_subset
otu_data <- OTU_Table_Fungi_Cluefield_Rhizosphere_subset

# metadata <- Phenotypes_2wpi_Cluefield
metadata <- Phenotypes_3wpi_Cluefield

otu_df <- otu_data %>% 
  select(-`#OTU ID`) %>% 
  group_by(taxonomy) %>% 
  summarise_all(list(sum)) %>% 
  gather(key = "SampleID", value = "Abund", -taxonomy) %>% 
  spread(key = "taxonomy", value = "Abund")

grep("taxonomy", colnames(otu_data))

#if the dataframe does not contain taxonomy info

otu_df <- otu_data %>% 
  # left_join(TaxInfo, by = "#OTU ID") %>% #For bacteria
  left_join(TaxInfoF, by = "#OTU ID") %>% #For fungi
  select(-`#OTU ID`) %>% 
  group_by(taxonomy) %>% 
  summarise_all(list(sum)) %>% 
  gather(key = "SampleID", value = "Abund", -taxonomy) %>% 
  spread(key = "taxonomy", value = "Abund")


dim(otu_df[complete.cases(otu_df),])
otu_df <- otu_df[complete.cases(otu_df),]
data.frame(Sample = otu_df$SampleID) %>% 
  separate(Sample, into = c("Plant", "Time", "Treat", "Replicate")) %>% View

xdata <- metadata %>% 
  filter(SampleID %in% otu_df$SampleID)
rownames(xdata) <- xdata$SampleID

otu_df_sub <- otu_df %>% 
  filter(SampleID %in% xdata$SampleID)
dim(otu_df_sub[complete.cases(otu_df_sub),])
rownames(otu_df_sub) <- otu_df_sub$SampleID

table(xdata[,3:5])

#Trimming the taxa abundance
library(gjam)
otu_df_subTr <- gjamTrimY(otu_df_sub[,-1], minObs = 3)
# otu_df_subTr <- gjamTrimY(otu_df_sub[,-1], minObs = 2) #Did this for the fungal endophytes
y <- otu_df_subTr$y
rownames(y) <- otu_df_sub$SampleID
dim(y)
tail(colnames(y))
xdata2 <- xdata
# xdata2$DMBQ <- scale(xdata2$DMBQ)
# xdata2$Vanillic_acid <- scale(xdata2$Vanillic_acid)
# barplot(apply(xdata2[,7:11],2, var)) #For week2
barplot(apply(xdata2[,7:9],2, var, na.rm = TRUE)) #For week3
hist(log(xdata2$Aerenchyma))
hist(log(xdata2$Suberin))
xdata2$Aerenchyma <- scale(xdata2$Aerenchyma)#For week3
xdata2$Suberin <- scale(xdata2$Suberin)#For week3
barplot(apply(xdata2[,7:9],2, var, na.rm = TRUE)) #For week3

barplot(apply(y,2, var))
# ydata <- cbind(xdata2[,7:11],y) #For week2
ydata <- cbind(xdata2[,7:9],y) #For week3

# typeNames    <- c(rep('CA',ncol(xdata2[,7:9])),
#                   rep('CON',ncol(xdata2[,10:11])),
#                   rep('CC',ncol(y)))   # composition count data
#For week 3
typeNames    <- c(rep('CON',ncol(xdata2[,7:9])),
                  #rep('CON',ncol(xdata2[,10:11])),
                  rep('CC',ncol(y)))   # composition count data


#Gjam model
# rl <- list(r = 5, N = 30)
# CCgroupsVar = c(rep(0,38),rep(1,ncol(y)),rep(2,ncol(yF)))
# CCgroupsVar = c(rep(0,36),rep(1,ncol(y)),rep(2,ncol(yF)))
(ng = 10000);(burnin =3000)
ml <- list(ng = ng, burnin =burnin, typeNames = typeNames)
# ml <- list(ng = 5000, burnin =1000, typeNames = typeNames[-42])
# xdata2 <- data.frame(xdata)
colnames(xdata2)
colnames(ydata)

xdata2$Soil <- factor(xdata2$Soil)
xdata2$Striga <- factor(xdata2$Striga)
summary(xdata2)

output1 <- gjam(~ Soil*Striga, data.frame(xdata2), data.frame(ydata), modelList = ml)
# output1 <- gjam(~ Soil*Striga, data.frame(xdata2), data.frame(ydata)[,-42], modelList = ml) #did this for fungal endophytes
View(output1$parameters$betaTable)
View(output2$parameters$betaTable)
View(output1$parameters$betaStandXWTable)
View(output1$parameters$betaStandXTable)

dim(output1$parameters$betaTable[output1$parameters$betaTable$sig95=="*",])
dim(output2$parameters$betaTable[output2$parameters$betaTable$sig95=="*",])

dim(outputBactEndo2week$parameters$betaTable[outputBactEndo2week$parameters$betaTable$sig95=="*",])

#Week2
# outputBactBK2week <- output1
# outputBactEndo2week <- output1
# outputBactRh2week <- output1

# outputFunBK2week <- output1
# outputFunEndo2week <- output1
# outputFunRh2week <- output1

#Week 3
# outputBactBK3week <- output1
# outputBactEndo3week <- output1
# outputBactEndoSand3week <- output1
# outputBactRh3week <- output1

# outputFunBK3week <- output1
# outputFunEndo3week <- output1
# outputFunEndoSand3week  <- output1
# outputFunRh3week <- output1

#Mega Model with everything

#Bacteria
otu_dataBkB <- OTU_Table_Bacteria_Cluefield_bulksoil_subset %>% 
  gather(key = "SampleID", value = "Bulk", -`#OTU ID`, -taxonomy)
otu_dataEnB <- OTU_Table_Bacteria_Cluefield_Endophytes_root_soil_subset %>% 
  gather(key = "SampleID", value = "Endo", -`#OTU ID`)
otu_dataRhB <- OTU_Table_Bacteria_Cluefield_Rhizosphere_subset %>% 
  gather(key = "SampleID", value = "Rhizo", -`#OTU ID`)

TabBacTot <- otu_dataBkB %>% 
  left_join(otu_dataEnB, by = c("SampleID", "#OTU ID")) %>% 
  left_join(otu_dataRhB, by = c("SampleID", "#OTU ID")) %>% 
  gather(key = "Compartment", value = "Abund", Bulk:Rhizo) %>% 
  group_by(SampleID, Compartment, taxonomy) %>% 
  summarise(AbundT = sum(Abund)) %>% 
  ungroup() %>% 
  spread(key = "taxonomy", value = "AbundT")
View(TabBacTot)
TabBacTot <- TabBacTot[complete.cases(TabBacTot),]
dim(TabBacTot)

metadataTot <- plyr::rbind.fill(Phenotypes_2wpi_Cluefield, Phenotypes_3wpi_Cluefield)
xdataTotBac2w <- TabBacTot[,1:2] %>% 
  left_join(Phenotypes_2wpi_Cluefield, by = "SampleID")
xdataTotBac2w <- xdataTotBac2w[complete.cases(xdataTotBac2w),]
View(xdataTotBac2w)
table(xdataTotBac2w[,c(2,5:6)])
xdataTotBac3w <- TabBacTot[,1:2] %>% 
  left_join(Phenotypes_3wpi_Cluefield, by = "SampleID")
xdataTotBac3w <- xdataTotBac3w[complete.cases(xdataTotBac3w),]
View(xdataTotBac3w)

#Fungi
otu_dataBkF <- OTU_Table_Fungi_Cluefield_bulksoil_subset
otu_dataEnF <- OTU_Table_Fungi_Cluefield_Endophytes_root_soil_subset
otu_dataRhF <- OTU_Table_Fungi_Cluefield_Rhizosphere_subset

TabFunTot <- otu_dataBkF %>% 
  left_join(otu_dataEnF, by = c("SampleID", "#OTU ID")) %>% 
  left_join(otu_dataRhF, by = c("SampleID", "#OTU ID")) %>% 
  gather(key = "Compartment", value = "Abund", Bulk:Rhizo) %>% 
  group_by(SampleID, Compartment, taxonomy) %>% 
  summarise(AbundT = sum(Abund)) %>% 
  ungroup() %>% 
  spread(key = "taxonomy", value = "AbundT")
View(TabFunTot)
TabFunTot <- TabFunTot[complete.cases(TabFunTot),]
dim(TabFunTot)
xdataTotBac <- TabFunTot[,1:2] %>% 
  left_join(Phenotypes_2wpi_Cluefield, by = "SampleID")
