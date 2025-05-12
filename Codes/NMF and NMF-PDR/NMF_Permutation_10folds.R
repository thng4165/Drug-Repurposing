

rm(list = ls())

setwd("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation")


#### I.Dataset
## Load data
load("rare_disease_interact.Rdata")
load("rare_disease_sim.RData")
load("rare_drug_sim.RData")

Asso_matrix = interact
disease_sim = rare_disease_sim
drug_sim = rare_drug_sim


# ## hsdn-MechDB
# load("hsdn_MechDB_dd_association.RData")
# load("hsdn_MechDB_drug_sim.RData")
# load("hsdn_MechDB_disease_sim_GIP.RData")
# Asso_matrix = hsdn_MechDB_dd_association
# disease_sim = hsdn_MechDB_disease_sim_GIP
# drug_sim = hsdn_MechDB_drug_sim




# # Cdata-MechDB
# load("Cdata_asso.RData")
# load("Cdata_drug_sim.RData")
# load("Cdata_disease_sim.RData")
# Asso_matrix = Cdata_asso
# disease_sim = Cdata_disease_sim
# drug_sim = Cdata_drug_sim
# 


## LAGCN dataset

# load("LAGCNdata_disease_sim.RData")
# load("LAGCNdata_drug_sim.RData")
# load("LAGCNdata_asso.RData")
# Asso_matrix = LAGCNdata_asso
# disease_sim = LAGCNdata_disease_sim
# drug_sim = LAGCNdata_drug_sim


## LRSSL dataset
# load("LRSSLdata_disease_sim.RData")
# load("LRSSLdata_drug_sim.RData")
# load("LRSSLdata_asso.RData")
# Asso_matrix = LRSSLdata_asso
# disease_sim = LRSSLdata_disease_sim
# drug_sim = LRSSLdata_drug_sim



## Ydataset
# load("Ydata_asso.RData")
# load("Ydata_drug_sim.RData")
# load("Ydata_disease_sim.RData")
# Asso_matrix = Ydata_asso
# disease_sim = Ydata_disease_sim
# drug_sim = Ydata_drug_sim

## DNdataset
# load("DNdata_asso.RData")
# load("DNdata_drug_sim.RData")
# load("DNdata_disease_sim.RData")
# Asso_matrix = DNdata_asso
# disease_sim = DNdata_disease_sim
# drug_sim = DNdata_drug_sim


#### II. functions

source("drugRepurposing_functions.R")


#### III. Data Preparation: CV for Association matrix and construct adjacency matrix

nFold = 10
#set.seed(2024)
## CV Association matrix
pos = c(Asso_matrix) > 0 
zer = c(Asso_matrix) == 0
length(pos) 
length(zer) 
pos_index = which(pos) 
zer_index = which(zer)
length(pos_index) 
length(zer_index) 


nr = nrow(Asso_matrix)  # Number of rows in the original matrix
nc = ncol(Asso_matrix)  # Number of columns in the original matrix



## Create nSeed seeds
nSeed = 1
# #list_seed = sample(1000:9999, nSeed, replace = FALSE)
# 
list_seed = 2024#c(4958, 8027, 6265, 9782, 1468, 5490, 7813, 8056, 5020, 9702, 2306, 1939, 4720, 3796, 9688, 8646, 8965,
#4985, 4270, 4915, 9148, 7998, 6512, 2167, 7835, 4008, 8873, 8864, 9799, 1483, 8853, 5772, 1764, 3988, 2303, 4773,
#4960, 1526, 5293, 2747, 6001, 9060, 7333, 1191, 5017, 7379, 2348, 3444, 6732, 7477)


AUC_null = rep(0, times=nSeed)
AUPR_null = rep(0, times=nSeed)



#for (iseed in 1 : nSeed){
iseed = 1

print(iseed)
  myseed = list_seed[iseed]
  

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
    fold[is.na(fold)] = 0 ## Initialize 0 for missing values of Association matrix for each fold
    
    temp = drug_sim %*% fold %*% disease_sim
    temp = fold * temp
    temp = as.matrix(temp)
    temp = temp/max(temp)
    CV_Adj_Input[[fold_index]] = temp
    
  }
  
  
  
  ###### IV. NMF Models
  results_list_element = list()
  set.seed(myseed)
  for( fold_index in 1:nFold){
    cat("\n processing the fold th ",fold_index)
    
    fold = Mising_Folds[[fold_index]] ## fold of Association matrix
    B0 = CV_Adj_Input[[fold_index]]  ## Adjacency matrix corresponding to fold
    x=B0[!is.na(fold)]
    x=x[x>0]
    theta=min(x)/(length(x)*100)
    B0[is.na(fold)]=theta
    B=B0
  
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
      train_set = fold
    )
    
    
    ############## V. Permutation NULL distribution
    o=result$test_approx
    ### make null distribution
    N=100
    k_ran=rank_fold
    useBoot=FALSE
    dat0=NULL
    set.seed(myseed)
    myseed2=sample(10000,1)
    
    cat("\n shuffling the adj values of train data :")
    for (id in 1:N){
      if (id %% 10 == 0) cat(" ",id)
      B_random = CV_Adj_Input[[fold_index]]
      pick=B_random >0 & !is.na(fold) #shuffle only ones with association=1
      B_random[pick] = sample(B_random[pick],replace = useBoot)
      B_random[is.na(fold)]=theta

      nmfB = NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
      B_update_random = nmfB$W %*% nmfB$H
      h0=B_update_random[is.na(fold)]
      dat0=rbind(dat0,h0)
    }
    shdat0=dat0 ##keep results of shuffling adj value
    result$shdat0=shdat0
    
    myseed2=sample(10000,1)
    dat0=NULL
    cat("\n sampling by rows :")
    for (id in 1:N){
      if (id %% 10 == 0) cat(" ",id)
      B_random =  CV_Adj_Input[[fold_index]]
      B_random = apply(B_random, 1, function(x){ #sampling by row
        testId=which(is.na(x))
        if (length(testId) > 0){
          x1=x[testId]
          x1[]= theta
          x2=x[-testId]
          x2=sample(x2, replace = useBoot)
          x[testId]=x1
          x[-testId]=x2
        }
        else{
          x=sample(x, replace = useBoot)
        }
        return(x)
      })

      B_random=t(B_random)

      nmfB = NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
      B_update_random = nmfB$W %*% nmfB$H
      h0=B_update_random[is.na(fold)]
      dat0=rbind(dat0,h0)
    }
    rdat0=dat0 ##keep results with sampling by rows
    result$rdat0=rdat0


    dat0=NULL
    myseed2=sample(10000,1)
    cat("\n sampling by columns :")
    for (id in 1:N){
      if (id %% 10 == 0) cat(" ",id)
      B_random =  CV_Adj_Input[[fold_index]]
      #B_random[is.na(fold)]=NA
      B_random = apply(B_random, 2, function(x){ #sampling by columns
        testId = which(is.na(x))
        if (length(testId) > 0){
          x1=x[testId]
          x1[]= theta
          x2=x[-testId]
          x2=sample(x2)
          x[testId]=x1
          x[-testId]=x2
        }
        else{
          x=sample(x)
        }
        return(x)
      })

      nmfB = NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2)
      B_update_random = nmfB$W %*% nmfB$H
      h0=B_update_random[is.na(fold)]
      dat0=rbind(dat0,h0)
    }
    cdat0=dat0 ##keep results with sampling by columns
    result$cdat0=cdat0

    
    dat0=NULL
    myseed2=sample(10000,1)
    cat("\n sampling by both rows and columns :")
    for (id in 1:N){
      if (id %% 10 == 0) cat(" ",id)
      B_random = c( CV_Adj_Input[[fold_index]])
      pick = is.na(B_random)
      B_random[!pick]=sample(B_random[!pick])
      B_random = matrix(B_random,nrow=nrow(fold),ncol=ncol(fold))
      B_random[is.na(fold)] = theta
      
      nmfB = NMF_ELee(V = B_random,  k = k_ran, itnum=100, eps=1e-9, epsDiff=1e-2, seed = myseed2) #standard
      B_update_random = nmfB$W %*% nmfB$H
      h0=B_update_random[is.na(fold)]
      dat0=rbind(dat0,h0)
    }
    bdat0=dat0 ##keep results with sampling by rows
    result$bdat0=bdat0


    result$o=o

    result$test_approx=o
    
    # Add the result to the list
    results_list_element[[fold_index]] = result
  }
  
  
  ####### VI. Predicted Matrix and Model's Performance
 
  ### Predicted matrix
  pred_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
  pr_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
  for (i in 1 : nFold){
    index = which(is.na(results_list_element[[i]]$train_set))
    pred_CV[index] = results_list_element[[i]]$test_approx
    
    result=results_list_element[[i]]

    o=result$o #if using standard NMF
    trueObs=Asso_matrix[index]
      
      
      
      shdat0=result$shdat0
      #dat0=rbind(result$cdat0)
      tstat_r=pr_r=result$o
      tstat_c=pr_c=result$o
      tstat_b=pr_b=result$o
      
      for (j in 1:length(o)){
        mytest=ks.test(shdat0[,j],result$rdat0[,j])
        pr_r[j]=mytest$p.value
        tstat_r[j]=mytest$statistic

        mytest=ks.test(shdat0[,j],result$cdat0[,j])
        pr_c[j]=mytest$p.value
        tstat_c[j]=mytest$statistic

        mytest=ks.test(shdat0[,j],result$bdat0[,j])
        pr_b[j]=mytest$p.value
        tstat_b[j]=mytest$statistic
      }
      
      
      o=tstat_r+tstat_c+tstat_b
      

       pred_CV[index] = o

  }
  

  ############### AUC - AUPR #######################
  
  inputObs_matrix=Asso_matrix
  prediction_matrix=pred_CV
  
  table(Observed=c(Asso_matrix)>0,Predicted=c(pred_CV)>0)

  
  #sort inputObs_matrix by column using the decreasing order by column of prediction_matrix
  res=sort_matrix(prediction_matrix,inputObs_matrix)
  
  sorted_inputObs_matrix=res$y_sorted
  sorted_score_matrix=res$score_sorted
  sort_index=res$sort_index
  
  #  sorted_inputObs_matrix=apply(sorted_inputObs_matrix,2,function(x) sample(x)) #make a random prediction, AUC should be 0.5
  
  tpr_list = NULL
  fpr_list = NULL
  recall_list = NULL
  precision_list = NULL
  accuracy_list = NULL
  F1_list = NULL
  
  #now compute performance metrics for top k=cutoff rows (which means top k predicted values for individual diseases or drugs in the columns) against the remaining rows
  for (cutoff in 1:nrow(sorted_inputObs_matrix)){
    P_matrix = sorted_inputObs_matrix[1:cutoff, ] #predicted Positives
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
  
  AUC=trapz(fpr_list,tpr_list)#AUC
  
  AUC
  # AUC_null[iseed] = AUC
  
  AUPR=trapz(recall_list,precision_list)#AUPR
  AUPR
  AUPR_null[iseed] = AUPR
  
#}
# save(AUC_null, file = "oMat_AUC_null_ttest_25runs.RData")
# save(AUPR_null, file = "oMat_AUPR_null_ttest_25runs.RData")

par(mfrow = c(1, 2))
plot(fpr_list,tpr_list, type="l", main=paste0("AUC=",round(AUC,3)), xlab="FPR (1-specificity)", ylab = "TPR (sensitivity)", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
abline(0,1,lty=2)
plot(recall_list,precision_list,type="l",main=paste0("AUPR=",round(AUPR,3)), xlab = "Recall", ylab="Precision", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
par(mfrow = c(1,1))



