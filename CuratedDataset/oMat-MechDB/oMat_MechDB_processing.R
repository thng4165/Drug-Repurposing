## Prepare data - WELL 

##------------------------------- Construct oMat-MechDB dataset -----------------
## Construct the oMat-MechDB dataset including three matrices: 
## 1. drug-disease association matrix, 2. drug-drug similarity matrix and 3. disease-disease similarity matrix
## such that diseases including in both oMat and MechDB
## Input databases: orphaNet (disease-symptom data), 
## In each step, datasets are saved to use later.


rm(list=ls())
setwd("C:/Trang/KIProjects/RareDiseases/oMatMechDB1/oMatMechDB/DataPreparation")


## --------------------disease-symptom data------------------------------------
## oMat provides the disease-symptom data.
load("Pre_datasets/orphaNet.RData")
head(oMat)
head(orphaNet_db)
length(oMat$Disease) ## 111765 instances
dim(oMat) ## 111765 x 5
#table(oMat$OrphaCode) ## the code and amount instances for each code
#table(oMat$Disease)   ## name of Disease and amount instances for each disease
unique_code = unique(oMat$OrphaCode) ## how many disease in data, not unique
num_unique_code = length(unique_code) ## 4240 diseases in oMat


##----------- Load orphaNet_Disease_mapping
## orphaNet_allIDs provides mesh id for disease in oMat data
## List of disease in oMat is subset of list of disease in orphaNet_Disease_mapping
load("Pre_datasets/orphanet_Disease_mapping.RData") # including "orphaNet_allIDs" and "orphaNet_synonyms"
head(orphaNet_allIDs)  # orphaNet_allIDs: inclding the id of diisease, each id presents as a feature of orphaNet_allIDs, meaning "id" is the same oMat$OrphaCode. For one id, there are including some information
head(orphaNet_synonyms) # one disease may be called by different names, but they should have the same id.
length(unique(orphaNet_synonyms))# 10839 with 
length(orphaNet_allIDs)
length(unique(orphaNet_allIDs)) ##7783

table(as.character(unique_code) %in% names(orphaNet_allIDs ))

##------------- prepare and extract disease_symptom_oMat

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

##------------ Check if each id has a unique name for oMat,
## names_per_id = TRUE, means that we have unique name for each Orphacode. 
## ===> We use orphaNet as a data of disease_symptoms
names_per_id <- data.frame(oMat) %>%
  group_by(OrphaCode) %>%
  summarise(unique_names = length(unique(Disease)))%>%
  pull(unique_names) %>%
  all()

## check missing values
sum(is.na(oMat$mesh)) ## 77028 indications with missing mesh in oMat

##---------- remove some NA mesh id, only NA in oMat
oMat_mesh = oMat$mesh[!is.na(oMat$mesh)] ## 1232 mesh id in oMat
length(unique(oMat_mesh))
## remove some empty mesh id ==> no blank space 
## length(unique(oMat_mesh)) = oMat_mesh[oMat_mesh!=""] 

##----------- Create Disease-symptom from orphaNet and orphaNet_Disease_Map
disease_symptom_oMat = oMat[oMat$mesh %in% oMat_mesh, ]
length(unique(disease_symptom_oMat$mesh)) ##1232 mesh id

length(unique(disease_symptom_oMat$Disease)) ## 1235 disease ==> synonym name
# save(disease_symptom_oMat, file = "Pre_datasets/disease_symptom_oMat.RData")


##-------------------------- DrugMechDB: drug-disease data ----------------------------------
## DrugMechDB provides drug-disease data
load("Pre_datasets/DrugMechDB.RData")
dim(DrugMechDB) ##4627 x 11


length(unique(DrugMechDB$drug_name)) ##1578
length(unique(DrugMechDB$drugbank)) ## 1558
length(unique(DrugMechDB$drug_mesh)) ## 1546

length(unique(DrugMechDB$disease_name)) ## 764
length(unique(DrugMechDB$disease_mesh)) ## 727


##------------ Prepare and extract drug-disease from MechDB data

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
# save(drug_disease_alls, file = "Pre_datasets/drug_disease_alls.RData")
head(drug_disease_alls)

##----------------------- Extract drugbank for drugs in drug-disease data (MechDB)
## chemInfo_DrugBank includes information and strucutre of drug. we use SMILES and INCHIKEY
load(file = "Pre_datasets/chemInfo_DrugBank.RData") 

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
# save(chemInfo_subset, file = "Pre_datasets/chemInfo_subset.RData")


## drug-disease-subset with information of drugbank
drug_disease_subset = drug_disease_alls[drug_disease_alls$drugbank %in% chemInfo_subset$DrugBankID, ]
any(unique(drug_disease_subset$drugbank) %in% chemInfo_subset$DrugBankID)
any( chemInfo_subset$DrugBankID %in% unique(drug_disease_subset$drugbank))
dim(drug_disease_subset)
# save(drug_disease_subset, file = "Pre_datasets/drug_disease_subset.RData")


## ----------------- Calculating Drug-Drug similarity Tanimoto for "subset" drugs in MechDB-----------------------

## load the dataset drug-disease
a = load("Pre_datasets/drug_disease_subset.RData")
load("Pre_datasets/chemInfo_subset.RData")
load("Pre_datasets/disease_symptom_oMat.RData")

## intersection of disease between disease-symptom-oMat, drug-disease-subset (small dataset)
## We did again in RareDiseaseDrug.R
intersect_disease = intersect(unique(disease_symptom_oMat$mesh), unique(drug_disease_subset$disease_mesh))
length(intersect_disease) # 89 diseases
## extract datasets with the overlap diseases between MechDB and oMat data
extract_disease_symptom = disease_symptom_oMat[disease_symptom_oMat$mesh %in% intersect_disease, ]
extract_drug_disease = drug_disease_subset[drug_disease_subset$disease_mesh %in% intersect_disease, ]
extract_chemInfo = chemInfo_subset[chemInfo_subset$DrugBankID %in% extract_drug_disease$drugbank, ]


## Drug-Drug similarity Tanimoto for drugs in MechDB after extracting drugbank, removing non-smiles, .....

# install.packages("BiocManager")
# BiocManager::install("ChemmineOB")
# install.packages("ChemmineR")
# install.packages("Cheminformatics")

library(ChemmineR)
library(Cheminformatics)


mysmiles = chemInfo_subset$SMILES
names(mysmiles) = chemInfo_subset$DrugBankID

sdfset = smiles2sdf(mysmiles) ### ==> 13 invalid SDFs detected ==> remove
valid = validSDF(sdfset)  ## check the valide entries
sdfset = sdfset[valid] ## only keep the valide entries and remove the invalide entries
apset = sdf2ap(sdfset) 
fpset = desc2fp(apset)
drug_sim = sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE)) ## drug-drug similarity
dim(drug_sim) ## 1319 x 1319 (remove 13 invalid SDFs strucutures)
drug_sim[1:5,1:5]

heatmap(drug_sim,
        col = colorRampPalette(c("white", "red"))(100),
        main = "Drug - Drug",
        xlab = "Drugs",
        ylab = "Drug",
        margins = c(5, 5))

# save(drug_sim, file = "Pre_datasets/drug_sim.RData") ## WELL DONE



##----------------- Extract sub-data of drug in MechDB again ---------------------------------
## Since 13 SMILES of 13 drugs have invalid SDF structures, and we remove them. it means that the chemInfo_subset
## is reduced 13 drugs. We extract again the drug_disease_subset then two these datasets have the same list of drugs

invalid_indices = which(!valid) ## indices of invalid term

## dataset after removing drug with invalid SMILES --> SDF
chemInfo_subset_s = chemInfo_subset[-invalid_indices, ] ## drugbank data after remove 
removed_DrugBank = chemInfo_subset$DrugBankID[invalid_indices] ## indicates of removed elements
drug_disease_subset_s = drug_disease_subset[!drug_disease_subset$drugbank %in% removed_DrugBank, ]
dim(drug_disease_subset_s) ## 4124 x 11 

#save(chemInfo_subset_s, file = "Pre_datasets/chemInfo_subset_s.RData") # after remove invalid validSDF
#save(drug_disease_subset_s, file="Pre_datasets/drug_disease_subset_s.RData") # after remove invalid validSDF
## WELL-DONE


##----------------------Calculating disease-disease similarity------------------------------------------
## Similarity of disease-disease by Gaussian Interaction Profile kernel, base on HPOId (symptoms)

##----------- Diseases from MechDB
load("Pre_datasets/drug_disease_subset_s.RData")
head(drug_disease_subset_s)
length(unique(drug_disease_subset_s$drugbank)) # 1332 drugs
length(unique(drug_disease_subset_s$disease_mesh)) # 674 diseases
length(unique(drug_disease_subset_s$disease_name)) # 705 diseases names

disease_MechDB = drug_disease_subset_s[, c("disease_name", "disease_mesh")]
disease_MechDB = unique(disease_MechDB)  #Remove duplicates


library(dplyr)
## remove the diseases have the same meshID, keep only the first occurrence
disease_MechDB = disease_MechDB %>% distinct(disease_mesh, .keep_all = TRUE)

# save(disease_MechDB, file = "Pre_datasets/disease_MechDB_s.RData")


##------ Diseases from oMat (disease - symptom)
load("Pre_datasets/disease_symptom_oMat.RData")
### Get unique diseases and HPO terms
unique_diseases <- unique(disease_symptom_oMat$mesh)
unique_hpo_terms <- unique(disease_symptom_oMat$HPOId)

## Create a matrix where rows represent diseases and columns represent HPO terms
intersection_matrix <- matrix(0, nrow = length(unique_diseases), ncol = length(unique_hpo_terms))

## Populate the interaction matrix based on the presence or absence of HPO terms for each disease
for (i in 1:length(unique_diseases)) {
  # i=1
  disease <- unique_diseases[i]
  hpo_terms <- unique(disease_symptom_oMat[disease_symptom_oMat$mesh == disease, "HPOId"])
  col_indices <- match(hpo_terms, unique_hpo_terms)
  intersection_matrix[i, col_indices] <- 1
}
dim(intersection_matrix) ## 1232 * 5535 = num_disease x num_HPOId


##------- Function to calculate Gaussian Interaction Profile kernel between two disease interaction profiles
gip_kernel <- function(profile1, profile2, gamma) {
  euclidean_distance <- sqrt(sum((profile1 - profile2)^2))
  similarity_value <- exp(-gamma * euclidean_distance)
  return(similarity_value)
}

#-------- Function to calculate disease-disease similarity matrix using GIP kernel
GIP_disease_matrix <- function(data_matrix, gamma) {
  num_diseases <- nrow(data_matrix)
  similarity_matrix <- matrix(NA, nrow = num_diseases, ncol = num_diseases)
  
  for (i in 1:num_diseases) {
    for (j in 1:num_diseases) {
      disease_profile1 <- data_matrix[i, ]
      disease_profile2 <- data_matrix[j, ]
      
      # Calculate GIP similarity
      similarity_matrix[i, j] <- gip_kernel(disease_profile1, disease_profile2, gamma)
    }
  }
  
  return(similarity_matrix)
}

##--------- Set the gamma parameter for the GIP kernel
hat_gamma = 1 ## Consider more for the choice of hat_gamma
num = nrow(intersection_matrix)
s = 0
for (i in 1:num) {
  disease_profile = intersection_matrix[i, ]
  s = s + sum(disease_profile^2)
}

gamma = hat_gamma/(1/num * s) ## 0.035

##------------ Calculating disease-disease similarity for diseases on oMat
## We may only calculate disease similarity for overlaped diseases between oMat and MechDB
disease_sim_GIP <- GIP_disease_matrix(intersection_matrix, gamma)
colnames(disease_sim_GIP) =  unique_diseases
rownames(disease_sim_GIP) =  unique_diseases
# Print the result
print(disease_sim_GIP[1:5, 1:5])

heatmap(disease_sim_GIP, 
        col = colorRampPalette(c("white", "red"))(100), # Color gradient
        main = "Disease Similarity Heatmap",
        xlab = "Diseases",
        ylab = "Diseases",
        margins = c(5, 5),
        labRow = unique_diseases,
        labCol = unique_diseases)
isSymmetric(disease_sim_GIP) ### TRUE
## The matrix is symmetric, but the heatplot looks nonsymmetric

# save(disease_sim_GIP, file = ("Pre_datasets/disease_sim_GIP.RData"))



##---------------------- Collecting final oMat-MechDB dataset ---------------------------------
##------------------ Drug disease matrix
## Picking diseases overlapping between oMat and MechDB after processing
## Picking drugs in MechDB with at least 1 association with overlapping diseases
## Getting drug-disease matrix for overlapping diseases between oMat and MechDB


load("Pre_datasets/chemInfo_subset_s.RData")
load("Pre_datasets/drug_disease_subset_s.RData")
load("Pre_datasets/disease_symptom_oMat.RData")

rare_intersection = intersect(disease_symptom_oMat$mesh, drug_disease_subset_s$disease_mesh)
drug_rare_disease = drug_disease_subset_s[drug_disease_subset_s$disease_mesh %in% rare_intersection, ]

## Matrix indicating the interactions between drugs and diseases in subset drug_rare_disease (same disease in two datasets)
interact = table(drug_rare_disease$drugbank, drug_rare_disease$disease_mesh)
interact = as.matrix(interact)

# save(interact, file = "Pre_datasets/rare_disease_interact.RData")

## -------------- Disease-disease similarity
load("Pre_datasets/disease_sim_GIP.RData") # disease-disease similarity of all diseases in oMat
dim(disease_sim_GIP) ## 1232 * 1232 Rare disease

## Picking disease-disease similarity corresponding to the overlapping diseases between oMat and MechDB
rare_disease_sim = disease_sim_GIP[rare_intersection,rare_intersection] ## 89 x 89
#save(rare_disease_sim, file = "Pre_datasets/rare_disease_sim.RData")

## ----------------- Drug-drug similarity
load("Pre_datasets/drug_sim.RData")
dim(drug_sim) ## 1319 * 1319

## Picking drug-drug similarity corresponding to drugs which have at least 1 association with overlapping diseases
rare_drug_sim = drug_sim[unique(drug_rare_disease$drugbank), unique(drug_rare_disease$drugbank)] # 148 * 148
#save(rare_drug_sim, file = "Pre_datasets/rare_drug_sim.RData")


# 
# ## check the order of the name of disease, drug. 
# any(unique(rownames(rare_drug_sim)) == unique(rownames(interact))) ##TRUE
# any(unique(rownames(rare_disease_sim)) == unique(colnames(interact))) ## TRUE

