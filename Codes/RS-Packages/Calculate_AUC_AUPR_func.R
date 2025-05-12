######################## Calculate AUC - AUPR ##############################

sort_matrix = function(score_matrix,interact_matrix){
  #For each disease in columns, sort the score values (predicted) and change the values of the interaction matrix accordingly. In downstream step, go down to each row, compute performance metrics for top 1, 2, 3, .. to summarize the final results.
  
  sort_index=score_matrix
  for ( i in 1:ncol(sort_index)){
    x=-score_matrix[,i]
    #x=score_matrix[,i]
    o=order(x, decreasing = FALSE)
    sort_index[,i]=o
  }
  
  score_sorted = score_matrix
  y_sorted = interact_matrix
  for (i in 1:ncol(interact_matrix)){
    score_sorted[,i] = score_matrix[,i][sort_index[,i]]
    y_sorted[,i] = interact_matrix[,i][sort_index[,i]]
  }
  rownames(score_sorted)=NULL #row names is not meaningful, so remove them to avoid potential misunderstanding
  rownames(y_sorted)=NULL
  
  res=list(y_sorted=y_sorted,score_sorted=score_sorted,sort_index=sort_index)
  return(res)
  
}



Compute_AUC_AUPR = function(inputObs_matrix,prediction_matrix){
  
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
  
  ##TRang add
  cut_list = NULL
  TP_list = NULL
  
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
    ## TRang Add
    cut_list = c(cut_list, cutoff)
    TP_list = c(TP_list, TP/88)
  }
  
  res_curve = list(tpr_list=tpr_list, fpr_list=fpr_list, recall_list=recall_list, precision_list=precision_list)
  return(res_curve)
}



AUC_AUPR = function(fpr_list, tpr_list, recall_list, precision_list){
  
  ### now compute AUC and AUPR
  library(caTools)
  
  AUC=trapz(fpr_list,tpr_list)#AUC
  AUC
  
  AUPR=trapz(recall_list,precision_list)#AUPR
  AUPR
  
  
  par(mfrow = c(1, 2))
  plot(fpr_list,tpr_list, type="l", main=paste0("AUC=",round(AUC,3)), xlab="FPR (1-specificity)", ylab = "TPR (sensitivity)", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2, xlim=c(0,1),ylim=c(0,1))
  abline(0,1,lty=2)
  plot(recall_list,precision_list,type="l",main=paste0("AUPR=",round(AUPR,3)), xlab = "Recall", ylab="Precision", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
  #plot(cut_list, TP_list, type="l", 
  #    main= "TP curve",
  #   xlab = "Top index",
  #  ylab = "True positive",
  # cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
  par(mfrow = c(1, 1))
  return (results = c(AUC, AUPR))
}
