### Here we calculate the disease-disease similarity matrix based on symptom

setwd("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom")

###############################################################################
## Here we calculate the disease-disease similarity of drug-disease dataset
## diseases in drug-disease dataset include name, mesh_id of disease

load("data/drug_disease_subset_s.RData")

head(drug_disease_subset_s)
length(unique(drug_disease_subset_s$drugbank)) # 1319 drugs
length(unique(drug_disease_subset_s$disease_mesh)) # 674 diseases
length(unique(drug_disease_subset_s$disease_name)) # 705 diseases names

### Extract the set of diseases from drug-disease data to calculate diseases distance
disease_MechDB = drug_disease_subset_s[, c("disease_name", "disease_mesh")]
disease_MechDB = unique(disease_MechDB)  #Remove duplicates


library(dplyr)
## remove the diseases have the same meshID, keep only the first occurrence
disease_MechDB = disease_MechDB %>% distinct(disease_mesh, .keep_all = TRUE)

save(disease_MechDB, file = "Data/disease_MechDB_s.RData")



################################################################################3

#### Similarity disease-disease matrix for disease in oMat based on symptom

load("Data/disease_symptom_oMat.RData")
######################## Disease-disease similarity GIP ###############################
### Similarity of disease-disease by Gaussian Interaction Profile kernel, base on HPOId also
### Get unique diseases and HPO terms

unique_diseases <- unique(disease_symptom_oMat$mesh)
unique_hpo_terms <- unique(disease_symptom_oMat$HPOId)

## Create a matrix where rows represent diseases and columns represent HPO terms
interaction_matrix <- matrix(0, nrow = length(unique_diseases), ncol = length(unique_hpo_terms))

## Populate the interaction matrix based on the presence or absence of HPO terms for each disease
for (i in 1:length(unique_diseases)) {
  i=1
  disease <- unique_diseases[i]
  hpo_terms <- unique(disease_symptom_oMat[disease_symptom_oMat$mesh == disease, "HPOId"])
  col_indices <- match(hpo_terms, unique_hpo_terms)
  interaction_matrix[i, col_indices] <- 1
}
dim(interaction_matrix) ## 1232 * 5535 = num_disease x num_HPOId
# Function to calculate Gaussian Interaction Profile kernel between two disease interaction profiles
gip_kernel <- function(profile1, profile2, gamma) {
  euclidean_distance <- sqrt(sum((profile1 - profile2)^2))
  similarity_value <- exp(-gamma * euclidean_distance)
  return(similarity_value)
}

# Function to calculate disease-disease similarity matrix using GIP kernel
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

# Set the gamma parameter for the GIP kernel
hat_gamma = 1 ## Consider more for the choice of hat_gamma
num = nrow(interaction_matrix)
s = 0
for (i in 1:num) {
  disease_profile = interaction_matrix[i, ]
  s = s + sum(disease_profile^2)
}

gamma = hat_gamma/(1/num * s) ## 0.035

# Apply the function to your disease interaction matrix
disease_sim_GIP <- GIP_disease_matrix(interaction_matrix, gamma)
colnames(disease_sim_GIP) =  unique_diseases
rownames(disease_sim_GIP) =  unique_diseases
# Print the result
print(S_disease_GIP[1:5, 1:5])

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


save(disease_sim_GIP, file = ("disease_sim_GIP.RData"))


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library("ComplexHeatmap")

Heatmap(disease_sim_GIP, show_column_names =FALSE,
        show_row_names = FALSE)


###################################Similarity disease-disease (Jaccard) ########################
### can not run for that, too slow
### Similarity of disease-disease by Jaccard base on symptom HPOId
# Install and load required packages
install.packages("stringdist")
install.packages("data.table")
install.packages("matrixStats")

library(stringdist)
library(data.table)
library(matrixStats)


## Calculate Jaccard similarity between two sets of strings
jaccard <- function(set1, set2) {
  intersect_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  return(intersect_size / union_size)
}

# Function to calculate disease-disease similarity matrix
unique_diseases <- unique(disease_symptom_oMat$mesh)
disease_sim_Jac <- matrix(NA, nrow = length(unique_diseases), ncol = length(unique_diseases))
colnames(disease_sim_Jac) <- unique_diseases
rownames(disease_sim_Jac) <- unique_diseases

for (i in 1:length(unique_diseases)) {
  for (j in 1:length(unique_diseases)) {
    disease1 <- unique_diseases[i]
    disease2 <- unique_diseases[j]
    
    hpo_set1 <- disease_symptom_oMat[disease_symptom_oMat$mesh == disease1, "HPOId"]
    hpo_set2 <- disease_symptom_oMat[disease_symptom_oMat$mesh == disease2, "HPOId"]
    
    # print(paste("Disease 1:", disease1, "HPO Set 1:", toString(hpo_set1)))
    # print(paste("Disease 2:", disease2, "HPO Set 2:", toString(hpo_set2)))
    
    
    disease_sim_Jac[i, j] <- jaccard(hpo_set1, hpo_set2)
  }
}

## save(S_disease_Jac,file="S_disease_Jac.RData")
## print similarity matrix
disease_sim_Jac[1:5, 1:5]
dim(disease_sim_Jac) ## 106 x 106. 

## Set the row and column names for better interpretation
rownames(disease_sim_Jac) <- colnames(disease_sim_Jac) <- unique_diseases

# Plot heatmap
heatmap(disease_sim_Jac,
        col = colorRampPalette(c("white", "red"))(100),
        main = "Disease Jaccard Similarity Heatmap",
        xlab = "Diseases",
        ylab = "Diseases",
        margins = c(5, 5),
        labRow = unique_diseases,
        labCol = unique_diseases)
