rm(list = ls())
setwd("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation")

myseed = 2024
set.seed(myseed)

#### I. Load data

load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_asso.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_drug_sim.RData")
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_disease_sim.RData")
Asso_matrix = t(HDVDdata_asso)
disease_sim = HDVDdata_disease_sim
drug_sim = HDVDdata_drug_sim

#load("Results/NMF_HDVDdata_process.RData")

load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/data_process/NMFPer_HDVD_process_kstest_25.RData")



#### II. Functions
source("drugRepurposing_functions.R")
nFold = 10
AUC = AUPR = NULL
for (iseed in 1:length(NMFLee_process)){
pred_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
Stat_CV= matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
pr_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
results_list_element = NMFLee_process[[iseed]]$NMF_process
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
    pr_r[j]=mytest$p.value
    
    # mytest=t.test(shdat0[,j],result$cdat0[,j])
    mytest=wilcox.test(shdat0[,j],result$cdat0[,j])
    # mytest=ks.test(shdat0[,j],result$cdat0[,j])
    tstat_c[j]=mytest$statistic
    pr_c[j]=mytest$p.value
    
    # mytest=t.test(shdat0[,j],result$bdat0[,j])
    mytest=wilcox.test(shdat0[,j],result$bdat0[,j])
    # mytest=ks.test(shdat0[,j],result$bdat0[,j])
    tstat_b[j]=mytest$statistic
    pr_b[j]=mytest$p.value
  }
  sumStats=tstat_r+tstat_c+tstat_b
  #sumStats = tstat_b
  fisherStats=-log10(pr_r)+-log10(pr_c)+-log10(pr_b)
  
  Stat_CV[index] = sumStats
  pr_CV[index] = fisherStats
}




inputObs_matrix=(Asso_matrix)

prediction_matrix = (Stat_CV)



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


### now compute AUC and AUPR
library(caTools)

NMF_AUC=trapz(fpr_list,tpr_list)#AUC
NMF_AUC


NMF_AUPR=trapz(recall_list,precision_list)#AUPR
NMF_AUPR

AUC[iseed] = NMF_AUC
AUPR[iseed] = NMF_AUPR

}


save(AUC, file = "NMFPer_AUC_HDVD_Wiltest_25.RData")
save(AUPR, file = "NMFPer_AUPR_HDVD_Wiltest_25.RData")

### List of AUC and AUPR
## NMF_Standard


## 1. HDVD - OK
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_HDVDdata_AUC_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_HDVDdata_AUPR_25.RData")

## 2. LAGCN - OK
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_LAGCN_25runs.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_LAGCN_25runs.RData")

## 3. Fdata - OK
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_Fdata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_Fdata_25.RData")

## 4. Cdata: ok
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_Cdata_25runs.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_Cdata_25runs.RData")

## 5. LRSSL - - OK
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_LRSSLdata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_LRSSLdata_25.RData")

## 6. Ydata
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_Ydata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_Ydata_25.RData")


## 7. SCMFDDL
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_SCMFDDLdata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_SCMFDDLdata_25.RData")

## 8. iDrug need to calculate AUC and AUPR again for iDrug data
## rerun
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_iDrug_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_iDrug_25.RData")

## 9. DNdata
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_DNdata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_DNdata_25.RData")


## 10. TLHGBI
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_TLHGBIdata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_TLHGBIdata_25.RData")

## 11. oMat - OK
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/oMat_AUC_100.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/oMat_AUPR_100.RData")

## 12. hsdn
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUC_hsdn_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_Standard/NMFS_AUPR_hsdn_25.RData")


## iDrug 






##### NMF Permutation t.test
## 4. Cdata
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/NMFPer_Cdata_ttest_AUC_100.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/NMFPer_Cdata_ttest_AUPR_100.RData")

load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/NMFPer_AUC_ttest_Cdata_25_rerun.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/NMFPer_AUPR_ttest_Cdata_25_rerun.RData")

## 1. HDVD
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/HDVD/NMFPer_HDVD_process_ttest_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/HDVD/NMFPer_AUC_HDVD_ttest_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/HDVD/NMFPer_AUPR_HDVD_ttest_25.RData")


## Fdata --ok

load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Fdata/NMFPer_AUC_ttest_Fdata_8.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Fdata/NMFPer_AUC_ttest_Fdata_8_16.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Fdata/NMFPer_AUC_ttest_Fdata_17_24.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Fdata/NMFPer_AUPR_ttest_Fdata_8.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Fdata/NMFPer_AUPR_ttest_Fdata_8_16.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Fdata/NMFPer_AUPR_ttest_Fdata_17_24.RData")


## Ydata --ok

load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Ydata/NMFPer_AUC_ttest_Ydata_8.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Ydata/NMFPer_AUC_ttest_Ydata_8_16.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Ydata/NMFPer_AUC_ttest_Ydata_17_24.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Ydata/NMFPer_AUPR_ttest_Ydata_8.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Ydata/NMFPer_AUPR_ttest_Ydata_8_16.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/Ydata/NMFPer_AUPR_ttest_Ydata_17_24.RData")


### LRSSL
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/LRSSL/NMFPer_ttest_LRSSL_AUC_AUPR_8_0.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/LRSSL/NMFPer_ttest_LRSSL_AUC_AUPR_8_1.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/LRSSL/NMFPer_ttest_LRSSL_AUC_AUPR_8_2.RData")

## SCMFDDL
# NMFLee_process[[4]]$AUPR_stat
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/SCMFDDL/NMFPer_ttest_SCMFDDL_AUC_AUPR_8_0.RData")



### DNdata
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/DNdata/NMFPer_AUC_ttest_DNdata_10.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/DNdata/NMFPer_AUPR_ttest_DNdata_10.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_ttest/DNdata/NMFPer_DNdata_AUC_AUPR_10.RData")

##### NMF Permutation ks.test
## 1. HDVD
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/HDVD/NMFPer_AUC_kstest_HDVD_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/HDVD/NMFPer_AUPR_kstest_HDVD_25.RData")

## 4. Cdata - dont choose kstest for Cdata
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/Cdata/NMFPer_AUC_kstest_Cdata_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/Cdata/NMFPer_AUPR_kstest_Cdata_25.RData")

## hsdn
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/hsdn/NMFPer_kstest_hsdn_AUC_AUPR_8_0.RData")


## idrug - results strange
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/iDrug/NMFPer_iDrug_kstest_AUC_AUPR_2_0.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_kstest/iDrug/NMFPer_iDrug_kstest_AUC_AUPR_8_6.RData")



##### NMF permutation wilcox test
## 2. LAGCN - OK 
#AUC
#AUPR
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/LAGCN/NMFPer_AUC_wiltest_LAGCN_8.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/LAGCN/NMFPer_AUC_wiltest_LAGCN_16.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/LAGCN/NMFPer_AUC_wiltest_LAGCN_17_24.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/LAGCN/NMFPer_AUPR_wiltest_LAGCN_8.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/LAGCN/NMFPer_AUPR_wiltest_LAGCN_16.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/LAGCN/NMFPer_AUPR_wiltest_LAGCN_17_24.RData")


### HDVD
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/HDVD/NMFPer_AUC_HDVD_Wiltest_25.RData")
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_25runs/NMF_per_wiltest/HDVD/NMFPer_AUPR_HDVD_Wiltest_25.RData")



#### HGIMC 

## HDVD
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_HDVD_AUC_25runs.csv")
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_HDVD_AUPR_25runs.csv")
## Cdata
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_Cdata_AUC_25runs.csv")
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_Cdata_AUPR_25runs.csv")
## Fdata
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_Fdata_AUC_25runs.csv")
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_Fdata_AUPR_25runs.csv")

## hsdn
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_hsdn_AUC_25runs.csv")
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_hsdn_AUPR_25runs.csv")

## LAGCN
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_LAGCN_AUC_25runs.csv")
read.csv("C:/Trang/KIProjects/ComprehensionDR/HGIMC/25runs/save/HGIMC_LAGCN_AUPR_25runs.csv")
