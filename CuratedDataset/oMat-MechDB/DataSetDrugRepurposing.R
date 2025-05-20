## Prepare data - WELL 

rm(list=ls())
setwd("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/oMatMechDB/DataPreparation/")

## This R script load the datasets which are used in project
## We also check the information of dataset.

#######################Load orphaNet#################################
## oMat provides the disease-symptom data.
load("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/orphaNet.RData")
head(oMat)
head(orphaNet_db)
length(oMat$Disease) ## 111765 instances
dim(oMat) ## 111765 x 5
#table(oMat$OrphaCode) ## the code and amount instances for each code
#table(oMat$Disease)   ## name of Disease and amount instances for each disease
unique_code <- unique(oMat$OrphaCode) ## how many disease in data, not unique
num_unique_code <- length(unique_code) ## 4240 diseases in oMat


########################### orphaNet_Disease_mapping.RData #############################
## orphaNet_allIDs provides mesh id for disease in oMat data
## List of disease in oMat is subset of list of disease in orphaNet_Disease_mapping
load("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/orphanet_Disease_mapping.RData") # including "orphaNet_allIDs" and "orphaNet_synonyms"
head(orphaNet_allIDs)  # orphaNet_allIDs: inclding the id of diisease, each id presents as a feature of orphaNet_allIDs, meaning "id" is the same oMat$OrphaCode. For one id, there are including some information
head(orphaNet_synonyms) # one disease may be called by different names, but they should have the same id.
length(unique(orphaNet_synonyms))# 10839 with 
length(orphaNet_allIDs)
length(unique(orphaNet_allIDs)) ##7783

table(as.character(unique_code) %in% names(orphaNet_allIDs ))

###################### prepare and extract disease_symptom_oMat #################

library(dplyr)
## Remove duplicate observations based on all columns
oMat = oMat %>% distinct()


## extract MeSH id from orphaNet_allIDs based on OrphaCode (same OrphaCode in both data), add mesh id in to oMat data. 
oMat$mesh = NA
#count = 0
id_extract = as.character(unique(oMat$OrphaCode))
for (i in 1 : length(id_extract)){
  id  = id_extract[i]
  t = which(names(orphaNet_allIDs) == id)
  x = orphaNet_allIDs[[t]]
  y = which(x$DB == "MeSH")
  if (length(y)>0){
    p = which(oMat$OrphaCode == id)
    oMat$mesh[p] = paste(x$ID[y],collapse = " " )
  }
#  count = count + length(y)
}

## check list of disease in oMat is subset of list disease of orparNet_allIDs by checking orphaCode.
id_extract_orpha = as.character(unique(names(orphaNet_allIDs)))
orphacode_id = intersect(id_extract,id_extract)
any(orphacode_id == id_extract) ## TRUE

## unlist mesh id for disease in oMat, ==> MAY WE DONT NEED HERE
## mesh_idalls = strsplit(oMat$mesh, " ")
## mesh_idalls = unlist(mesh_idalls)  
## oMat$mesh = mesh_idalls ## not work

## Check if each id has a unique name for oMat, names_per_id = TRUE, means that we have unique name for each Orphacode. 
## ===> We use orphaNet as a data of disease_symptoms
names_per_id <- data.frame(oMat) %>%
  group_by(OrphaCode) %>%
  summarise(unique_names = length(unique(Disease)))%>%
  pull(unique_names) %>%
  all()

## check missing values
sum(is.na(oMat$mesh)) ## 77028 indications with missing mesh in oMat

## remove some NA mesh id, only NA in oMat
oMat_mesh = oMat$mesh[!is.na(oMat$mesh)] ## 1232 mesh id in oMat
length(unique(oMat_mesh))
## remove some empty mesh id ==> no blank space 
## length(unique(oMat_mesh)) = oMat_mesh[oMat_mesh!=""] 

## Create Disease-symptom from orphaNet and orphaNet_Disease_Map
disease_symptom_oMat = oMat[oMat$mesh %in% oMat_mesh, ]
length(unique(disease_symptom_oMat$mesh)) ##1232 mesh id

length(unique(disease_symptom_oMat$Disease)) ## 1235 disease ==> synonym name
save(disease_symptom_oMat, file = "data/disease_symptom_oMat.RData")


############################ DrugMechDB ####################################
## DrugMechDB provides drug-disease data
load("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/DrugMechDB.RData")
dim(DrugMechDB) ##4627 x 11


length(unique(DrugMechDB$drug_name)) ##1578
length(unique(DrugMechDB$drugbank)) ## 1558
length(unique(DrugMechDB$drug_mesh)) ## 1546

length(unique(DrugMechDB$disease_name)) ## 764
length(unique(DrugMechDB$disease_mesh)) ## 727


############################# Prepare and extract drug-disease#######################

## Remove duplicate observations based on all columns
DrugMechDB = DrugMechDB %>% distinct()

## rewrite mesh id in DrugMechDB data - same format in oMat to find the keyword between two data
disease_mesh_DB = sapply(DrugMechDB$disease_mesh, function(x){strsplit(x, ":")[[1]][2]}) ## mesh id disease from drug-disease data
disease_mesh_DB = unlist(disease_mesh_DB)
length(unique(disease_mesh_DB)) ## 727 mesh diseases in drug-diseases data
length(unique(DrugMechDB$disease_name)) ## 764 diseases name ==> synonym name or missing mesh id
DrugMechDB$disease_mesh = disease_mesh_DB
length(unique(DrugMechDB$disease_mesh))
length(unique(DrugMechDB$drug_mesh)) ## 1546 drugs
drug_disease_alls = DrugMechDB

## remove some NA mesh id ==> no NA values for entire data
sum(is.na(drug_disease_alls))
## remove some empty mesh id ==> 61 drug_mesh are empty + other columns have few empty (blank space) => 75
length(which(drug_disease_alls == "")) ## ==> no 
drug_disease_alls = drug_disease_alls[rowSums(drug_disease_alls == "") == 0,] ## remove rows with at least 1 empty element
dim(drug_disease_alls) ## 4552
length(unique(drug_disease_alls$drug_mesh)) # 1545 drugs
save(drug_disease_alls, file = "data/drug_disease_alls.RData")
head(drug_disease_alls)

#################################### Extract drugbank for drug_disease_alls #######################

## chemInfo_DrugBank includes information and strucutre of drug. we use SMILES and INCHIKEY
load(file = "D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/chemInfo_DrugBank.RData") 

## intersection of drugbank between chemInfo_DrugBank and drug_disease_alls
subset_drugbank = intersect(chemInfo_DrugBank$DrugBankID, drug_disease_alls$drugbank) #1496
## extract the subset of chemInfo with subset drugbank
chemInfo_subset = chemInfo_DrugBank[chemInfo_DrugBank$DrugBankID %in% subset_drugbank, ]## 1496 x x25
head(chemInfo_subset)

## SMILES information
mysmiles = as.character(chemInfo_subset$SMILES)
names(mysmiles) = chemInfo_subset$DrugBankID
mysmiles=mysmiles[!is.na(mysmiles)] #remove some without SMILEs information
mysmiles=mysmiles[mysmiles!=""]#remove some without SMILEs information
length(mysmiles) ## 1332 after remove some rows without SMILES

## extract subset of cheInfo after remove some without SMILES
chemInfo_subset = chemInfo_subset[chemInfo_subset$DrugBankID %in% names(mysmiles), ]
dim(chemInfo_subset) ## 1332 x 24
save(chemInfo_subset, file = "data/chemInfo_subset.RData")

######################### drug-disease-subset datasets #############################
drug_disease_subset = drug_disease_alls[drug_disease_alls$drugbank %in% chemInfo_subset$DrugBankID, ]
any(unique(drug_disease_subset$drugbank) %in% chemInfo_subset$DrugBankID)
any( chemInfo_subset$DrugBankID %in% unique(drug_disease_subset$drugbank))
dim(drug_disease_subset)
save(drug_disease_subset, file = "data/drug_disease_subset.RData")



