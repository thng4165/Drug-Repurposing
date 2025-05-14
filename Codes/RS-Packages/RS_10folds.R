
rm(list = ls())
setwd("RSmodels")

myseed = 123456
set.seed(myseed)

library(recommenderlab)

#### I. Load data
## oMat data
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/rare_disease_interact.Rdata")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/rare_disease_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/rare_drug_sim.RData")
Asso_matrix = interact
Asso_matrix = as.matrix(Asso_matrix) ## table
disease_sim = rare_disease_sim
drug_sim = rare_drug_sim


## HDVD data
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_disease_sim.RData")
Asso_matrix = t(HDVDdata_asso)
disease_sim = HDVDdata_disease_sim
drug_sim = HDVDdata_drug_sim

## LAGCN
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_disease_sim.RData")
Asso_matrix = (LAGCNdata_asso)
disease_sim = LAGCNdata_disease_sim
drug_sim = LAGCNdata_drug_sim

## Fdata
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Fdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Fdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Fdata_disease_sim.RData")
Asso_matrix = Fdata_asso
disease_sim = Fdata_disease_sim
drug_sim = Fdata_drug_sim


## Cdata
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Cdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Cdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Cdata_disease_sim.RData")
Asso_matrix = Cdata_asso
disease_sim = Cdata_disease_sim
drug_sim = Cdata_drug_sim

## Ydata
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Ydata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Ydata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Ydata_disease_sim.RData")
Asso_matrix = Ydata_asso
disease_sim = Ydata_disease_sim
drug_sim = Ydata_drug_sim



## hsdn-MechDB
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/hsdn_MechDB_dd_association.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/hsdn_MechDB_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/hsdn_MechDB_disease_sim_GIP.RData")
Asso_matrix = hsdn_MechDB_dd_association
disease_sim = hsdn_MechDB_disease_sim_GIP
drug_sim = hsdn_MechDB_drug_sim

## LRSSL
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LRSSLdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LRSSLdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LRSSLdata_disease_sim.RData")
Asso_matrix = (LRSSLdata_asso)
disease_sim = LRSSLdata_disease_sim
drug_sim = LRSSLdata_drug_sim


#### II. Functions
source("drugRepurposing_functions.R")


#### III. Data Preparation: CV for Association matrix and construct adjacency matrix
nFold = 10
## CV Association matrix
pos = c(Asso_matrix) > 0 
zer = c(Asso_matrix) == 0
pos_index = which(pos) 
zer_index = which(zer)


nr = nrow(Asso_matrix)  # Number of rows (drugs) in the original matrix
nc = ncol(Asso_matrix)  # Number of columns (diseases) in the original matrix
# fold the data
pos_fold = makeFolds(k = nFold, n = length(pos_index)) ## positive
zer_fold = makeFolds(k = nFold, n = length(zer_index)) ## zeros


# Combine positive and negative folds alternately (ten folds)
ele_fold = vector("list", length = length(pos_fold))
for (i in seq_along(pos_fold)) {
  ele_fold[[i]] = c(pos_index[pos_fold[[i]]], zer_index[zer_fold[[i]]])
}


## convert element in fold by NA in original matrix
## Folds include 0, positive, NA 
## Folds- of Association matrix with missing values
Mising_Folds = lapply(seq_along(ele_fold), function(fold_index){
  inde = ele_fold[[fold_index]]
  fold_MM = Asso_matrix ## initial for the convert matrix
  p_list = list() ## list of position
  
  for (i in 1:length(inde)) {
    #print(paste("Loop iteration:", i, "Value of i in inde:", inde[i]))
    p = convert_element(nr, nc, k = inde[i])
    p_list[[length(p_list) + 1]] = p
    fold_MM[p[1], p[2]] = NA
  }
  return (as.matrix(fold_MM))
})


### adj matrix of each fold inlcudes the missing values from the cv Asso matrix
Ad_cv = list()
Ad_mat = list()
for (i in 1: nFold){
  Asso = Mising_Folds[[i]]
  Ad_cv[[i]] = matrix(0, nrow=nr+nc, ncol=nr+nc) # adjacency matrix for each fold
  row.names(Ad_cv[[i]]) = c(rownames(drug_sim),rownames(disease_sim))
  colnames(Ad_cv[[i]]) = c(colnames(drug_sim), colnames(disease_sim))
  Ad_cv[[i]][1:nr, 1:nr] = as.matrix(drug_sim)
  Ad_cv[[i]][1:nr,(nr+1):(nr+nc)] = Asso
  Ad_cv[[i]][(nr+1):(nr+nc), 1:nr] = t(Asso)
  Ad_cv[[i]][(nr+1):(nr+nc), (nr+1):(nr+nc)] = as.matrix(disease_sim)
  Ad_mat[[i]] = as.matrix(Ad_cv[[i]])
  Ad_mat[[i]] = as(Ad_mat[[i]], "realRatingMatrix")
  # Ad_mat[[i]] = normalize(Ad_mat[[i]])

}
### Prediction
pred_matrix = list()
pred_asso = list()
for (i in 1: nFold){
    ## IBCF package
 # reco = Recommender(Ad_mat[[i]], "IBCF",param=list(normalize = NULL, method="Cosine"))
    ## LIBMF
    reco = Recommender(Ad_mat[[i]], "LIBMF")  # )# param=list(normalize = NULL, method="Cosine"))
  
  predictions = predict(reco, Ad_mat[[i]], type="ratings")
  pred_matrix[[i]] = as.matrix(predictions@data)
  pred_asso[[i]] = pred_matrix[[i]][1:nr,(nr+1):(nr+nc)]
}



## check if the predictions of two parts associations equal to each other
for (i in 1:nFold){
  a = any(pred_matrix[[i]][1:nr,(nr+1):(nr+nc)]==t(pred_matrix[[i]][(nr+1):(nr+nc), 1:nr]))
  print(a)
}

## Final predicted score matrix
pred_score =  matrix(NA,nrow = nr, ncol = nc)
for (i in 1:nFold) {
  index = which(is.na(Mising_Folds[[i]]))
  pred_score[index] = pred_asso[[i]][index]
}


####################################################
## AUC and AUPR
source("Calculate_AUC_AUPR_func.R")

inputObs_matrix = Asso_matrix
prediction_matrix = pred_score
res = Compute_AUC_AUPR(inputObs_matrix,prediction_matrix)
library(caTools)

AUC_AUPR_result = AUC_AUPR(res$fpr_list,res$tpr_list,res$recall_list,res$precision_list)

