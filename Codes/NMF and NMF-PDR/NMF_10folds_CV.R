rm(list = ls())

setwd("C:/Trang/KIProjects/ComprehensionDR/RerunMatlabCode/NMF_NMFPermutation")


start.time <- Sys.time()


myseed = 2024
set.seed(myseed)

#### I. Load data
# load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/rare_disease_interact.Rdata")
# load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/rare_disease_sim.RData")
# load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/rare_drug_sim.RData")
# Asso_matrix = (interact)
# disease_sim = rare_disease_sim
# drug_sim = rare_drug_sim

load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_disease_sim.RData")
Asso_matrix = t(HDVDdata_asso)
disease_sim = HDVDdata_disease_sim
drug_sim = HDVDdata_drug_sim

load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_disease_sim.RData")
Asso_matrix = (LAGCNdata_asso)
disease_sim = LAGCNdata_disease_sim
drug_sim = LAGCNdata_drug_sim


# load("Fdata_asso.RData")
# load("Fdata_drug_sim.RData")
# load("Fdata_disease_sim.RData")
# Asso_matrix = Fdata_asso
# disease_sim = Fdata_disease_sim
# drug_sim = Fdata_drug_sim



load("Cdata_asso.RData")
load("Cdata_drug_sim.RData")
load("Cdata_disease_sim.RData")
Asso_matrix = Cdata_asso
disease_sim = Cdata_disease_sim
drug_sim = Cdata_drug_sim



## hsdn-MechDB
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/hsdn_MechDB_dd_association.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/hsdn_MechDB_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/hsdn_MechDB_disease_sim_GIP.RData")
Asso_matrix = hsdn_MechDB_dd_association
disease_sim = hsdn_MechDB_disease_sim_GIP
drug_sim = hsdn_MechDB_drug_sim

load("LRSSLdata_asso.RData")
load("LRSSLdata_drug_sim.RData")
load("LRSSLdata_disease_sim.RData")
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


nr = nrow(Asso_matrix)  # Number of rows in the original matrix
nc = ncol(Asso_matrix)  # Number of columns in the original matrix
# fold the data
pos_fold = makeFolds(k = nFold, n = length(pos_index), seedVal = myseed) ## positive
zer_fold = makeFolds(k = nFold, n = length(zer_index), seedVal = myseed) ## zeros


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
  return (fold_MM)
})

## Construct the adjacency matrix corresponding 10 folds of Association matrix
CV_Adj_Input = list()

for (fold_index in 1:nFold){
  fold = Mising_Folds[[fold_index]]
  fold[is.na(fold)] = 0
  # fold[is.na(fold)] = BNNR[is.na(fold)]
  temp = drug_sim %*% fold %*% disease_sim
  temp = fold * temp
  temp = as.matrix(temp)
  temp = temp/max(temp)
  
  CV_Adj_Input[[fold_index]] = temp
}


###### IV. NMF Models
results_list_element = list()

for (fold_index in 1: nFold){
# fold_index = 1  
  cat("\n processing the fold th ",fold_index)
  
  fold = Mising_Folds[[fold_index]] ## fold of Association matrix
  
  B0 = CV_Adj_Input[[fold_index]]  ## Adjacency matrix corresponding to fold
  x=B0[!is.na(fold)]
  x=x[x>0]
  theta=min(x)/(length(x)*100)
  B0[is.na(fold)]=theta
  B=B0
  
  
  set.seed(myseed)
  ## rank of B - since we perform nmf with B
  rank_fold = find_optimal_rank(B0, 0.90) ## contain 95%
  
  nmfB = NMF_ELee(V = B, k = rank_fold, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed) #standard
  
  ## NMF approximation
  B_update = nmfB$W %*% nmfB$H
  
  
  
  ## Store the results in a list
  result = list(
    fold_number = fold_index,
    train_approx = B_update,
    test_approx = B_update[is.na(fold)],
    M2 = B,
    train_set = fold
  )
  results_list_element[[fold_index]] = result
  
  
  
  
  ############## V. Permutation NULL distribution
  o=result$test_approx
  ### make null distribution
  N=100
  k_ran=rank_fold # use the same rank of the standard NMF
  
  
  ###
  useBoot=FALSE
  dat0=NULL
  set.seed(myseed)
  myseed2=sample(10000,1)
  
  shdat0=rdat0=cdat0=bdat0=NULL
  
  shdat1=rdat1=cdat1=bdat1=NULL
  for (id in 1:N){
    #cat("\n shuffling the adj values of train data :")
    if (id %% 10 == 0) cat(" ",id)
    B_random = CV_Adj_Input[[fold_index]]
    pick=B_random >0 & !is.na(fold) #shuffle only ones with association=1
    
    B_random[pick] = sample(B_random[pick],replace = useBoot)
    B_random[is.na(fold)]= theta 
    nmfB <- NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
    B_update_random = nmfB$W %*% nmfB$H
    h0=B_update_random[is.na(fold)]
    shdat0=rbind(shdat0,h0)
    #shdat1=rbind(shdat1,B_update_random[!is.na(fold)])
    
    
    #cat("\n sampling by rows :")
    B_random =  CV_Adj_Input[[fold_index]]
    B_random[is.na(fold)]= theta
    B_random = apply(B_random, 1, function(x) sample(x))
    B_random=t(B_random) #need to be transposed
    
    nmfB <- NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
    B_update_random = nmfB$W %*% nmfB$H
    h0=B_update_random[is.na(fold)]
    rdat0=rbind(rdat0,h0)
    #rdat1=rbind(rdat1,B_update_random[!is.na(fold)])
    
    #cat("\n sampling by columns :")
    B_random =  CV_Adj_Input[[fold_index]]
    # B_random[is.na(fold)]= theta
    B_random = apply(B_random, 2, function(x) sample(x))
    nmfB <- NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
    B_update_random = nmfB$W %*% nmfB$H
    h0=B_update_random[is.na(fold)]
    cdat0=rbind(cdat0,h0)
    #cdat1=rbind(cdat1,B_update_random[!is.na(fold)])
    
    #cat("\n sampling by both rows and columns :")
    B_random=CV_Adj_Input[[fold_index]]
    # B_random[is.na(fold)]= theta
    B_random=sample(c(B_random))
    B_random=matrix(B_random,nrow=nrow(fold),ncol=ncol(fold))
    nmfB <- NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
    B_update_random = nmfB$W %*% nmfB$H
    h0=B_update_random[is.na(fold)]
    bdat0=rbind(bdat0,h0)
    #bdat1=rbind(bdat1,B_update_random[!is.na(fold)])
  }
  
  result$shdat0=shdat0
  result$rdat0=rdat0
  result$cdat0=cdat0
  result$bdat0=bdat0
  
  
  result$o=o
  results_list_element[[fold_index]] = result
}

pred_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
Stat_CV= matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
#pr_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])

for (i in 1 : nFold){
  
  cat("\n Processing fold ",i)
  index = which(is.na(results_list_element[[i]]$train_set))
  pred_CV[index] = results_list_element[[i]]$test_approx
  
  result=results_list_element[[i]]
  
  o=result$test_approx #if using standard NMF
  trueObs=Asso_matrix[index]
  
  shdat0=result$shdat0
  tstat_r=pr_r=result$test_approx
  tstat_c=pr_c=result$test_approx
  tstat_b=pr_b=result$test_approx
  
  
  for (j in 1:length(o)){
     # mytest=t.test(shdat0[,j],result$rdat0[,j])
    mytest=wilcox.test(shdat0[,j],result$rdat0[,j])
     # mytest=ks.test(shdat0[,j],result$rdat0[,j])
    tstat_r[j]=mytest$statistic
   # pr_r[j]=mytest$p.value
    
     # mytest=t.test(shdat0[,j],result$cdat0[,j])
    mytest=wilcox.test(shdat0[,j],result$cdat0[,j])
     # mytest=ks.test(shdat0[,j],result$cdat0[,j])
    tstat_c[j]=mytest$statistic
   # pr_c[j]=mytest$p.value
    
     # mytest=t.test(shdat0[,j],result$bdat0[,j])
    mytest=wilcox.test(shdat0[,j],result$bdat0[,j])
     # mytest=ks.test(shdat0[,j],result$bdat0[,j])
    tstat_b[j]=mytest$statistic
   # pr_b[j]=mytest$p.value
  }
  sumStats=tstat_r+tstat_c+tstat_b
 # fisherStats=-log10(pr_r)+-log10(pr_c)+-log10(pr_b)
  
  Stat_CV[index] = sumStats
 # pr_CV[index] = fisherStats
}

## HDVD update
# save(results_list_element, file = "C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/ResultsNMFPer/NMF_HDVD_process.RData")
# save(results_list_element, file = "C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/ResultsNMFPer/NMF_LAGCN_process_rerun.RData")


## Fdata with 
# save(results_list_element, file = "NMF_Fdata_process.RData")
# save(Stat_CV, file = "NMF_Per_Fdata_prediction.RData")
# save(pred_CV, file = "NMF_standard_Fdata_prediction.RData")


## Cdata
#save(results_list_element, file = "NMF_Cdata_process.RData")
# save(Stat_CV, file = "NMF_Per_Cdata_ttest_prediction.RData")
#save(Stat_CV, file = "NMF_Per_Cdata_kstest_prediction.RData")
#save(pred_CV, file = "NMF_standard_Cdata_prediction.RData")
#save(pr_CV, file = ("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest/NMF_Per_Pr_kstest_Cdata.RData"))




## oMat with ks.test 
# save(results_list_element, file = "NMF_oMatdata_process.RData")
# save(Stat_CV, file = "NMF_Per_oMatdata_prediction.RData")
# save(pred_CV, file = "NMF_standard_oMatdata_prediction.RData")

# ## HDVD with t.test provides better results than ks.test, and B_random[is.na(fold)] = 0 befor sampling obs cases
# save(results_list_element, file = "NMF_HDVDdata_process.RData")
# save(Stat_CV, file = "NMF_Per_HDVDdata_prediction.RData")
# save(pred_CV, file = "NMF_standard_HDVDdata_prediction.RData")


# ## LAGCN 
# save(results_list_element, file = "NMF_LAGCNdata_process.RData")
# save(Stat_CV, file = "NMFPer_LAGCNdata_Stat_kstest_prediction.RData")
# save(pr_CV, file = "NMFPer_LAGCNdata_pr_kstest_prediction.RData")
# save(pred_CV, file = "NMF_standard_LAGCNdata_prediction.RData")

## hsdndata with 
# save(results_list_element, file = "NMF_hsdndata_process.RData")
# save(Stat_CV, file = "NMF_Per_hsdndata_prediction.RData")
# save(pred_CV, file = "NMF_standard_hsdndata_prediction.RData")


## hsdndata with 
# save(results_list_element, file = "NMF_hsdndata_process.RData")
# save(Stat_CV, file = "NMF_Per_hsdndata_prediction.RData")
# save(pred_CV, file = "NMF_standard_hsdndata_prediction.RData")

inputObs_matrix=(Asso_matrix)

prediction_matrix = pred_CV
prediction_matrix = Stat_CV
prediction_matrix = pr_CV


#sort inputObs_matrix by column using the decreasing order by column of prediction_matrix
res=sort_matrix(prediction_matrix,inputObs_matrix) 

sorted_inputObs_matrix=res$y_sorted
sorted_score_matrix=res$score_sorted
sort_index=res$sort_index

tpr_list = NULL
fpr_list = NULL

recall_list = NULL
precision_list = NULL

accuracy_list = NULL
F1_list = NULL

#now compute performance metrics for top k=cutoff rows (which means top k predicted values for individual diseases or drugs in the columns) against the remaining rows
for (cutoff in 1:nrow(sorted_inputObs_matrix)){
  #for (cutoff in 1:length(sorted_inputObs_matrix)){
  P_matrix = sorted_inputObs_matrix[1:cutoff,] #predicted Positives
  N_matrix = sorted_inputObs_matrix[-c(1:cutoff),] #predicted Negatives
  
  TP = sum(P_matrix > 0)
  FP = sum(P_matrix == 0)
  TN = sum(N_matrix == 0)
  FN = sum(N_matrix > 0)
  tpr = TP / (TP + FN)
  fpr = FP / (FP + TN)
  recall = TP / (TP + FN)
  precision = TP / (TP + FP)
  accuracy = (TN + TP) / (TN + TP + FN + FP)
  F1 = (2 * TP) / (2*TP + FP + FN)
  if ((2*TP + FP + FN)==0)  F1 = 0
  
  F1_list=c(F1_list,F1)
  accuracy_list=c(accuracy_list,accuracy)
  
  tpr_list = c(tpr_list,tpr)
  fpr_list = c(fpr_list,fpr)
  recall_list = c(recall_list, recall)
  precision_list = c(precision_list,precision)
}


# LAGCN = list()
# LAGCN$fpr_list = fpr_list
# LAGCN$tpr_list = tpr_list
# LAGCN$recall_list = recall_list
# LAGCN$precision_list = precision_list
# save(LAGCN, file = "LAGCN_matrix_wilttest_NMFPer.RData")
# 
### now compute AUC and AUPR
library(caTools)

NMF_AUC=trapz(fpr_list,tpr_list)#AUC
NMF_AUC

NMF_AUPR=trapz(recall_list,precision_list)#AUPR
NMF_AUPR


par(mfrow = c(1, 2))
plot(fpr_list,tpr_list, type="l", main=paste0("AUC=",round(NMF_AUC,3)), xlab="FPR (1-specificity)", ylab = "TPR (sensitivity)", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
abline(0,1,lty=2)
plot(recall_list,precision_list,type="l",main=paste0("AUPR=",round(NMF_AUPR,3)), xlab = "Recall", ylab="Precision", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
par(mfrow = c(1,1))


