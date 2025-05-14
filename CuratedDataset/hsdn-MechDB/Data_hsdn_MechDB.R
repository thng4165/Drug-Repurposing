### Disease-symptom dataset from hsdn and MechDB
### Constructing Disease similarity matrix with symptom terms (GIP)
### Constructing Drug similarity matrix 
### Constructing Drug-disease association matrix

rm(list=ls())

# setwd("DrugSymptom/Mesh_hsdn/DataPreparation")


####################################Disease-Symptom dataset##################
load("disease_MechDB_s.RData") ## subset from MechDB data
load("hsdn_withMeshId.RData") ## hsdn data



### Adding the mesh id of disease into hsdn_net
### Dataset includes intersection diseases of hsdn_net and drug-disease MechDB

hsdn_disease_symptom = hsdn_net

hsdn_disease_symptom$order = seq_len((nrow(hsdn_disease_symptom)))
hsdn_disease_symptom = merge(hsdn_disease_symptom, hsdn_disease[, c("MeSH Disease Term", "Mesh_id")], by.x ="MeSH Disease Term", by.y="MeSH Disease Term", all.x = TRUE)
hsdn_disease_symptom = hsdn_disease_symptom[order(hsdn_disease_symptom$order), ]

#k = which(is.na(hsdn_disease_symptom$Mesh_id))
#hsdn_disease_symptom$`MeSH Disease Term`[k] ## check disease without Mesh_id
### 110 missing values in Mesh_id, corresponding to 1 disease without Mesh_id.

### remove the missing value (remove disease without Mesh_id)
hsdn_disease_symptom = na.omit(hsdn_disease_symptom)
### ===> 4218 disease, 4216 Mesh_id, 322 terms of symptom with 147868 observations 
### ===> 3 diseases have the same Mesh_id (synonym name)
### We work with Mesh_id, have not considered synonym names here



### Intersection disease set between MechDB data and hsdn
int = intersect(hsdn_disease_symptom$Mesh_id, disease_MechDB$disease_mesh)
### ===> 616 opverlaping diseases

### sub-data of hsdn disease-symptom with opvelaping diseases 
hsdn_MechDB = hsdn_disease_symptom[hsdn_disease_symptom$Mesh_id %in% int]
hsdn_MechDB = unique(hsdn_MechDB)


### Check the synonym names
x=paste0(hsdn_MechDB$`MeSH Symptom Term`,"__",hsdn_MechDB$Mesh_id)
y=which(duplicated(x))
p=which(x %in% x[y])
hsdn_MechDB[p,] ## ==> "Stevens-Johnson Syndrome" and "Epidermal Necrolysis, Toxic" are synonym names with the same symptoms, then we may remove one of them, e.g Epidermal Necrolysis, Toxic.

hsdn_MechDB = hsdn_MechDB[hsdn_MechDB$`MeSH Disease Term` != "Epidermal Necrolysis, Toxic", ]
## ==> no synonym names now
## ==> 35306 observation with 616 diseases

save(hsdn_MechDB, file = "hsdn_MechDB.RData")


#################################### Disease Similarity Matrix GIP ###################

unique_diseases = unique(hsdn_MechDB$Mesh_id)
unique_hpo = unique(hsdn_MechDB$`MeSH Symptom Term`)

## Create a matrix where rows represent diseases and columns represent HPO terms
interaction_matrix = matrix(0, nrow = length(unique_diseases), ncol = length(unique_hpo))


## Populate the interaction matrix based on the presence or absence of HPO terms for each disease
for (i in 1:length(unique_diseases)) {
  disease = unique_diseases[i]
  hpo_terms = unique(hsdn_MechDB$`MeSH Symptom Term`[hsdn_MechDB$Mesh_id == disease])
  col_indices = match(hpo_terms, unique_hpo)
  interaction_matrix[i, col_indices] = 1
}


# Function to calculate Gaussian Interaction Profile kernel between two disease interaction profiles
gip_kernel = function(profile1, profile2, gamma) {
  euclidean_distance = sqrt(sum((profile1 - profile2)^2))
  similarity_value = exp(-gamma * euclidean_distance)
  return(similarity_value)
}

# Function to calculate disease-disease similarity matrix using GIP kernel
GIP_disease_matrix = function(data_matrix, gamma) {
  num_diseases = nrow(data_matrix)
  similarity_matrix = matrix(NA, nrow = num_diseases, ncol = num_diseases)
  
  for (i in 1:num_diseases) {
    for (j in 1:num_diseases) {
      disease_profile1 = data_matrix[i, ]
      disease_profile2 = data_matrix[j, ]
      # Calculate GIP similarity
      similarity_matrix[i, j] = gip_kernel(disease_profile1, disease_profile2, gamma)
    }
  }
  return(similarity_matrix)
}

# Set the gamma parameter for the GIP kernel
hat_gamma = 1 ## Consider more for the choice of hat_gamma
num = nrow(interaction_matrix)
s = 0
for (i in 1:num) {
  disease_profile = interaction_matrix[i, ]
  s = s + sum(disease_profile^2)
}

gamma = hat_gamma/(1/num * s) ## 0.017

# Calculate the disease - symptom interaction matrix
hsdn_MechDB_disease_sim_GIP = GIP_disease_matrix(interaction_matrix, gamma)
colnames(hsdn_MechDB_disease_sim_GIP) =  unique_diseases
rownames(hsdn_MechDB_disease_sim_GIP) =  unique_diseases


save(hsdn_MechDB_disease_sim_GIP, file = ("hsdn_MechDB_disease_sim_GIP.RData"))


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")

png(file = "hsdn_MechDB_disease_sim_GIP.png")

Heatmap(hsdn_MechDB_disease_sim_GIP, show_column_names =FALSE,
        show_row_names = FALSE)
dev.off()


###################Disease Similarity Jacard distance ######################

#install.packages("stringdist")
#install.packages("data.table")
#install.packages("matrixStats")

library(stringdist)
library(data.table)
library(matrixStats)


## Calculate Jaccard similarity between two sets of strings
jaccard = function(set1, set2) {
  intersect_size = length(intersect(set1, set2))
  union_size = length(union(set1, set2))
  return(intersect_size / union_size)
}

# Function to calculate disease-disease similarity matrix
unique_diseases = unique(hsdn_MechDB$Mesh_id)
hsdn_MechDB_disease_sim_Jac = matrix(NA, nrow = length(unique_diseases), ncol = length(unique_diseases))
colnames(hsdn_MechDB_disease_sim_Jac) = unique_diseases
rownames(hsdn_MechDB_disease_sim_Jac) = unique_diseases

for (i in 1:length(unique_diseases)) {
  for (j in 1:length(unique_diseases)) {
    disease1 = unique_diseases[i]
    disease2 = unique_diseases[j]
    hpo_set1 = hsdn_MechDB$`MeSH Symptom Term`[hsdn_MechDB$Mesh_id == disease1]
    hpo_set2 = hsdn_MechDB$`MeSH Symptom Term`[hsdn_MechDB$Mesh_id == disease2]
    hsdn_MechDB_disease_sim_Jac[i, j] = jaccard(hpo_set1, hpo_set2)
  }
}

save(hsdn_MechDB_disease_sim_Jac,file="hsdn_MechDB_disease_sim_Jac.RData")
 
library("ComplexHeatmap")
png(file = "hsdn_MechDB_disease_sim_Jac.png")
Heatmap(hsdn_MechDB_disease_sim_Jac, show_column_names =FALSE,
        show_row_names = FALSE)
dev.off()


############################## drug similarity matrix ########################

#### Drug similarity matrix for all drug in MechDB after processing (clearning, check Smiles structure, valids...)

load("chemInfo_subset_s.RData") ## 1319 drug from intersection MechBD and Drugbank with SMILES structures.


#install.packages("BiocManager")
#BiocManager::install("ChemmineOB")
#install.packages("ChemmineR")
##install.packages("Cheminformatics")
#library(Cheminformatics)

library(ChemmineR)


mysmiles = chemInfo_subset_s$SMILES
names(mysmiles) = chemInfo_subset_s$DrugBankID

sdfset = smiles2sdf(mysmiles) 
valid = validSDF(sdfset) 
sdfset = sdfset[valid] 
apset = sdf2ap(sdfset) 
fpset = desc2fp(apset)
drug_sim = sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE)) ## drug-drug similarity

save(drug_sim, file = "MechDB_drug_sim.RData") 

library("ComplexHeatmap")
png(file = "MechDB_drug_sim.png")
Heatmap(drug_sim, show_column_names =FALSE,
        show_row_names = FALSE)
dev.off()


### extract drug similarity based on the intersection between MechDB and hsdn
load("hsdn_MechDB.RData") 
load("drug_disease_subset_s.RData")

### extract drugs from MechDB based on intersection diseases between MEchDB and hsdn
hsdn_MechDB_drug = drug_disease_subset_s
hsdn_MechDB_drug = hsdn_MechDB_drug[hsdn_MechDB_drug$disease_mesh %in% hsdn_MechDB$Mesh_id,] ## intersection by mesh id of disease ==> corresponding to 1270/1319 drugs.
hsdn_MechDB_drug_sim = drug_sim[unique(hsdn_MechDB_drug$drugbank), unique(hsdn_MechDB_drug$drugbank)] # 1270*1270

save(hsdn_MechDB_drug_sim, file = "hsdn_MechDB_drug_sim.RData")


######################### Drug-Disease Association matrix#######################
### hsdn_MechDB_drug includes 3701 observations ==> 3701 associations ==> 0.47% positive values
### 
hsdn_MechDB_dd_association = table(hsdn_MechDB_drug$drugbank, hsdn_MechDB_drug$disease_mesh)
save(hsdn_MechDB_dd_association, file = "hsdn_MechDB_dd_association.RData")






