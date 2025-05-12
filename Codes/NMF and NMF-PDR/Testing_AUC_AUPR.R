
rm(list = ls())
setwd("C:/Trang/KIProjects/ComprehensionDR/NMF_Lee_check_model")

#### I.Dataset
## Load data
load("rare_disease_interact.Rdata")
load("rare_disease_sim.RData")
load("rare_drug_sim.RData")
Asso_matrix = interact
disease_sim = rare_disease_sim
drug_sim = rare_drug_sim
# 

# ## hsdn-MechDB
# load("hsdn_MechDB_dd_association.RData")
# load("hsdn_MechDB_drug_sim.RData")
# load("hsdn_MechDB_disease_sim_GIP.RData")
# Asso_matrix = hsdn_MechDB_dd_association
# disease_sim = hsdn_MechDB_disease_sim_GIP
# drug_sim = hsdn_MechDB_drug_sim




## Cdataset
load("Cdata_asso.RData")
load("Cdata_drug_sim.RData")
load("Cdata_disease_sim.RData")
Asso_matrix = Cdata_asso
disease_sim = Cdata_disease_sim
drug_sim = Cdata_drug_sim


## Fdataset- 593 * 313

# load("Fdata_disease_sim.RData")
# load("Fdata_drug_sim.RData")
# load("Fdata_asso.RData")
# Asso_matrix = Fdata_asso
# disease_sim = Fdata_disease_sim
# drug_sim = Fdata_drug_sim
# 


## Ydataset

# load("Ydata_disease_sim.RData")
# load("Ydata_drug_sim.RData")
# load("Ydata_asso.RData")
# Asso_matrix = Ydata_asso
# disease_sim = Ydata_disease_sim
# drug_sim = Ydata_drug_sim


## LRSSL dataset
# load("LRSSLdata_disease_sim.RData")
# load("LRSSLdata_drug_sim.RData")
# load("LRSSLdata_asso.RData")
# Asso_matrix = LRSSLdata_asso
# disease_sim = LRSSLdata_disease_sim
# drug_sim = LRSSLdata_drug_sim


## LAGCN dataset
# load("LAGCNdata_disease_sim.RData")
# load("LAGCNdata_drug_sim.RData")
# load("LAGCNdata_asso.RData")
# Asso_matrix = LAGCNdata_asso
# disease_sim = LAGCNdata_disease_sim
# drug_sim = LAGCNdata_drug_sim

## DNdataset
# load("DNdata_asso.RData")
# load("DNdata_drug_sim.RData")
# load("DNdata_disease_sim.RData")
# Asso_matrix = DNdata_asso
# disease_sim = DNdata_disease_sim
# drug_sim = DNdata_drug_sim


## CV 10-folds

## NMF standard: prediction_matrix of models

#load("pred_CV_NMF_oMat.RData") 
# load("pred_CV_NMF_hsdn.RData") 
#load("pred_CV_NMF_Cdata.RData")
#load("pred_CV_NMF_Fdata.RData")
#load("pred_CV_NMF_LRSSLdata.RData") 
#load("pred_CV_NMF_Ydata.Rdata") 
#load("pred_CV_NMF_NDdata.RData") 


### Test 1 in slide (NMF permutation using empirical p-value and with smapling by row)
#load("pred_CV_Null_dist_oMat.RData") # test 1
# #load("pred_CV_Null_dist_hsdn.RData") # test 1
# load("pred_CV_Null_dist_oMat.RData") # test 1


### Test 2 in slide (NMF permutation with ttest)
#load("oMatdata_10fold_propNull_ttest_13_5.RData") # predicted value without Na for row and column shuffle
load("Cdata_10fold_propNull_ttest_13_5.RData") # predicted value without Na for row and column shuffle
##load("hsdndata_10fold_propNull_ttest_13_5.RData") # predicted value withour Na for row and column shuffle
# load("Fdata_10fold_propNull_ttest_13_5.RData")# predicted value without Na for row and column shuffle
#load("LRSSLdata_10fold_propNull_ttest_15_5.RData")
#load("Ydata_10fold_propNull_ttest_15_5.RData")
#load("DNdata_10fold_propNull_ttest_15_5.RData")


# rm(list=ls())
# R = list()
# result2 = load("DNdata_10fold_withNullDistr_ttest_18_5_1_8.RData")
# for (i in 1:8){
#   R[[i]] = results_list_element[[i]]
# }
# result1 = load("DNdata_10fold_propNull_ttest_18_5_9_19.RData")
# R[[9]] = results_list_element[[9]]
# R[[10]] = results_list_element[[10]]
# 
# for (i in 1:10){
#   results_list_element[[i]] = R[[i]]
# }
####save(results_list_element, file = "DNdata_10fold_propNull_ttest_15_5.RData")





## new Null distribution with rank_null = 2
 #load("oMatdata_10fold_newnull_ttest_13_5.RData")

## LOOCV drugs



## LOOCV diseases
#a = load("oMatdata_propNull_ttest_LOOCV_disease.RData") ## with B_random[is.na(fold)] = NA



source("drugRepurposing_functions.R")

### Predicted matrix
pred_CV_SumVal = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
pred_CV_SumFisher = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
pred_CV_ComVal = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
pred_CV_ComFisher = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
pr_CV = matrix(NA, nrow = dim(Asso_matrix)[1], ncol = dim(Asso_matrix)[2])
nr = nrow(Asso_matrix)  # Number of rows in the original matrix
nc = ncol(Asso_matrix)  # Number of columns in the original matrix



nFold = 10
for (i in 1 : nFold){
  index = which(is.na(results_list_element[[i]]$train_set))
  pred_CV_ComVal[index] = results_list_element[[i]]$test_approx
  pred_CV_ComFisher[index] = results_list_element[[i]]$test_approx

  result=results_list_element[[i]]

  o=result$test_approx #if using standard NMF
  result$o = o
  trueObs=Asso_matrix[index]


  #
  shdat0 = result$shdat0
  dat0=rbind(result$cdat0)
  tstat_r=pr_r=result$o
  tstat_c=pr_c=result$o
  tstat_b=pr_b=result$o
  zer = result$o-result$o
  o_tst = result$o
  o_pr = result$o

   for (j in 1:length(o)){
    mytest=t.test(shdat0[,j],result$rdat0[,j])
    #mytest=wilcox.test(shdat0[,j],result$rdat0[,j])
    #mytest=ks.test(shdat0[,j],result$rdat0[,j])
    pr_r[j]=mytest$p.value
    tstat_r[j]=mytest$statistic

    mytest=t.test(shdat0[,j],result$cdat0[,j])
    #mytest=wilcox.test(shdat0[,j],result$cdat0[,j])
    #mytest=ks.test(shdat0[,j],result$cdat0[,j])
    pr_c[j]=mytest$p.value
    tstat_c[j]=mytest$statistic

    mytest=t.test(shdat0[,j],result$bdat0[,j])
    #mytest=wilcox.test(shdat0[,j],result$bdat0[,j])
    #mytest=ks.test(shdat0[,j],result$bdat0[,j])
    pr_b[j]=mytest$p.value
    tstat_b[j]=mytest$statistic

    mytest=t.test(3*shdat0[,j]-result$rdat0[,j]-result$cdat0[,j]-result$bdat0[,j],zer )
    #mytest=ks.test(3*shdat0[,j]-result$rdat0[,j]-result$cdat0[,j]-result$bdat0[,j],zer )
    o_tst[j] = mytest$statistic
    o_pr[j] = mytest$p.value
  }


  o=tstat_r+tstat_c+tstat_b
  o_Fisher=-2*(log(pr_r)+log(pr_c)+log(pr_b))
  #o = pr_r + pr_c + pr_b

    o1 = o_tst
    o2 = o_pr


   pred_CV_SumVal[index] = o1
   pred_CV_SumFisher[index] = o2 

  pred_CV_ComVal[index] = o
  pred_CV_ComFisher[index] = o_Fisher
  
}

# save(pred_CV_SumFisher, file = "Cdata_PredVal_SumFisher_ttest.RData")
# save(pred_CV_SumVal, file = "Cdata_PredVal_SumVal_ttest.RData")
# save(pred_CV_ComVal, file = "Cdata_PredVal_ComVal_ttest.RData")
# save(pred_CV_ComFisher, file = "Cdata_PredVal_ComFisher_ttest.RData")


load("Cdata_PredVal_SumFisher_ttest.RData")
load("Cdata_PredVal_SumVal_ttest.RData")
load("Cdata_PredVal_ComVal_ttest.RData")
load("Cdata_PredVal_ComFisher_ttest.RData")

load("Cdata_PredVal_SumFisher_kstest.RData")
load("Cdata_PredVal_SumVal_kstest.RData")
load("Cdata_PredVal_ComVal_kstest.RData")
load("Cdata_PredVal_ComFisher_kstest.RData")


# load("oMat_PredVal_ComVal_test2.RData")
true_values = Asso_matrix[index]
sumStat_Val = pred_CV_SumVal[index]
sumStat_Pr = pred_CV_SumFisher[index]
comTest_Val = pred_CV_ComVal[index]
comTest_Pr = pred_CV_ComFisher[index]





true_values_kstest = Asso_matrix[index]
sumStat_Val_kstest = pred_CV_SumVal[index]
sumStat_Pr_kstest = pred_CV_SumFisher[index]
comTest_Val_kstest = pred_CV_ComVal[index]
comTest_Pr_kstest = pred_CV_ComFisher[index]



plot(train_data, train_standard)
test_pred = matrix(0,nrow = nrow(train_standard), ncol = ncol(train_standard))

for (i in 1 : 10){
  index = which(is.na(results_list_element[[i]]$train_set))
  test_pred[index] = results_list_element[[i]]$test_approx
}

i = 1
result = results_list_element[[i]]
index = which(is.na(result$train_set))
train_standard = result$train_approx
test_standard = result$test_approx
train_data = result$train_set
train_data[is.na(train_data)] = 0

true_values_ttest = Asso_matrix[index]
sumStat_Val_ttest = pred_CV_SumVal[index]
sumStat_Pr_ttest = pred_CV_SumFisher[index]
comTest_Val_ttest = pred_CV_ComVal[index]
comTest_Pr_ttest = pred_CV_ComFisher[index]



trainPred=result$train_approx[-index]
trainGroup=Asso_matrix[-index]
prepdata1=data.frame(var1=trainPred, var2= trainStat, Type=trainGroup)



plot(prepdata1$var1, prepdata1$var2,  pch=trainPch,  xlim=c(0,max(stat1)), xlab="NMF prediction", ylab="statistic")
# plot(prepdata$var1, prepdata$var2,  pch=trainPch, xlab="NMF prediction", ylab="statistic")
# points(result$o[trueObs!=1],sumStats[trueObs!=1], pch=mypch+16, col="blue",cex=1.5, lwd=2)
# points(result$o[trueObs==1],sumStats[trueObs==1], pch=mypch+18, col="red",cex=1.5, lwd=2)
# legend("bottom",c("asso=1","asso=0"),pch=c(8,1))



plot(train_standard, test_pred, main= "Predicted values: NMF standard",
     ylab="Predicted test values",
     xlab = "NMF prediction")
points( train_standard[index][true_values_ttest == 0], test_pred[index][true_values_ttest == 0], pch = 1, col = "blue")
points( train_standard[index][true_values_ttest == 1], test_pred[index][true_values_ttest == 1], pch = 8, col = "red")
legend("bottomright", legend = c("True Value 0", "True Value 1"),
       col = c("blue", "red"), pch = c(1, 8))




a = load("oMat_BNNR_FullMatrix_Null.RData")

test_pred = matrix(0,nrow = nrow(train_standard), ncol = ncol(train_standard))

for (i in 1 : 10){
  index = which(is.na(results_list_element[[i]]$train_set))
  test_pred[index] = results_list_element[[i]]$test_approx
}
test_pred = t(test_pred)

i = 1


result = results_list_element[[i]]
index = which(is.na(result$train_set))
true_values = t(Asso_matrix)[index]

train_standard = result$train_approx
test_standard = result$test_approx
train_data = result$train_set
train_data[is.na(train_data)] = 0


plot(c((train_standard)), main= "BNNR prediction",
     ylab="BNNR prediction",
     xlab = "index")
points(index, ((train_standard)[index]), pch = 19, col = "blue")
points(which(true_values == 1),train_standard[which(true_values == 1)], pch = 24, col = "red")


plot(((train_standard[-index])), test_pred[-index], main= "BNNR prediction",
     ylab="BNNR prediction",
     xlab = "index")
points( c((train_standard)[index]), pch = 19, col = "blue")
points( c((train_standard)[true_values == 1]), pch = 24, col = "red")



legend("bottomright", legend = c("Asso = 0", "Asso = 1"),
       col = c("blue", "red"), pch = c(24, 19))




plot(t(test_pred)[index], main= "BNNR prediction: test dataset",
     ylab="BNNR prediction",
     xlab = "index")
points(which(true_values == 0), t(test_pred)[index][which(true_values == 0)], pch = 19, col = "blue")
points(which(true_values == 1), t(test_pred)[index][which(true_values == 1)], pch = 24, col = "red")

legend("topright", legend = c("Asso = 0", "Asso = 1"),
       col = c("blue", "red"), pch = c(19, 24))


plot(c(train_standard[-index]), main= "BNNR prediction: train dataset",
     ylab="BNNR prediction",
     xlab = "index")

points( which((true_values) == 0), c(test_standard)[which((true_values) == 0)], pch = 19, col = "blue")
points(which((true_values) == 1), c(test_standard)[which((true_values) == 1)], pch = 24, col = "red")

legend("topright", legend = c("Asso = 0", "Asso = 1"),
       col = c("blue", "red"), pch = c(19, 24))




# 
# plot(comTest_Val, sumStat_Val, main= "Predicted values (t-statistics): NMF permutation", 
#      ylab="t-statistics: Combining 3 tests", 
#      xlab = "t-statistics: testing for sum of observed and null data")
# points( comTest_Val[true_values == 0], sumStat_Val[true_values == 0], pch = 1, col = "blue")
# points( comTest_Val[true_values == 1], sumStat_Val[true_values == 1], pch = 8, col = "red")
# legend("bottomright", legend = c("True Value 0", "True Value 1"),
#        col = c("blue", "red"), pch = c(1, 8))
# 
# 
# plot(comTest_Pr, sumStat_Pr, main= "P-values: NMF permutation", 
#      ylab="P-values: Combining 3 tests", 
#      xlab = "P-Values: testing for sum of observed and null data")
# points( comTest_Pr[true_values == 0], sumStat_Pr[true_values == 0], pch = 1, col = "blue")
# points( comTest_Pr[true_values == 1], sumStat_Pr[true_values == 1], pch = 8, col = "red")
# legend("topright", legend = c("True Value 0", "True Value 1"),
#        col = c("blue", "red"), pch = c(1, 8))
# 
# 
# 
plot(result$train_approx[index], sumStat_Val_ttest, main= "NMF standard and NMF permutation",
     ylab="t-statistics of NMF permutation",
     xlab = "Predicted values of NMF standard")
points(result$train_approx[index][true_values == 0], sumStat_Val_ttest[true_values == 0], pch = 1, col = "blue")
points(result$train_approx[index][true_values == 1], sumStat_Val_ttest[true_values == 1], pch = 8, col = "red")
legend("bottomright", legend = c("True Value 0", "True Value 1"),
       col = c("blue", "red"), pch = c(1, 8))


plot(result$train_approx[index], sumStat_Val_kstest, main= "NMF standard and NMF permutation",
     ylab="ks-statistics of NMF permutation",
     xlab = "Predicted values of NMF standard")
points(result$train_approx[index][true_values == 0], sumStat_Val_kstest[true_values == 0], pch = 1, col = "blue")
points(result$train_approx[index][true_values == 1], sumStat_Val_kstest[true_values == 1], pch = 8, col = "red")
legend("bottomright", legend = c("True Value 0", "True Value 1"),
       col = c("blue", "red"), pch = c(1, 8))
# 



plot(sumStat_Val_ttest, sumStat_Val_kstest, main= "NMF standard and NMF permutation",
     xlab="t-statistics",
     ylab = "ks statistics")
points(sumStat_Val_ttest[true_values == 0], sumStat_Val_kstest[true_values == 0], pch = 1, col = "blue")
points(sumStat_Val_ttest[true_values == 1], sumStat_Val_kstest[true_values == 1], pch = 8, col = "red")
legend("bottomright", legend = c("True Value 0", "True Value 1"),
       col = c("blue", "red"), pch = c(1, 8))
# 
# 
# 
# plot(result$train_approx[index], sumStat_Pr, main= "NMF standard and NMF permutation", 
#      ylab="P-Values of NMF permutation", 
#      xlab = "Predicted values of NMF standard")
# points(result$train_approx[index][true_values == 0], sumStat_Pr[true_values == 0], pch = 1, col = "blue")
# points(result$train_approx[index][true_values == 1], sumStat_Pr[true_values == 1], pch = 8, col = "red")
# legend("topright", legend = c("True Value 0", "True Value 1"),
#        col = c("blue", "red"), pch = c(1, 8))
# 
# 
# plot(result$train_approx[index], comTest_Pr, main= "NMF standard and NMF permutation", 
#      #ylab="t-statistics of NMF permutation", 
#      ylab="p-values of NMF permutation",
#      xlab = "Predicted values of NMF standard")
# points(result$train_approx[index][true_values == 0], comTest_Pr[true_values == 0], pch = 1, col = "blue")
# points(result$train_approx[index][true_values == 1], comTest_Pr[true_values == 1], pch = 8, col = "red")
# legend("topright", legend = c("True Value 0", "True Value 1"),
#        col = c("blue", "red"), pch = c(1, 8))
# 
# 
# 
# 
# plot(result$train_approx[index], comTest_Val, main= "NMF standard and NMF permutation", 
#      ylab="t-statistics of NMF permutation", 
#      xlab = "Predicted values of NMF standard")
# points(result$train_approx[index][true_values == 0], comTest_Val[true_values == 0], pch = 1, col = "blue")
# points(result$train_approx[index][true_values == 1], comTest_Val[true_values == 1], pch = 8, col = "red")
# legend("topright", legend = c("True Value 0", "True Value 1"),
#        col = c("blue", "red"), pch = c(1, 8))
# 
# 
# 
 i = 1
 
 j = 10
 
 result=results_list_element[[j]]
 plot(result$o)
 
 plot(results_list_element[[10]]$test_approx)
 
 
 result$train_set
 pick = which(is.na(result$train_set))
 
 c = Asso_matrix[pick]
d = which(c<1)    

e = which.max(result$test_approx[d])

 
a = result$rdat0[,i]

par(mfrow = c(1, 1))



pick_st = which.max((tstat_r))

i = 822
plot(result$shdat0[,i], main = "Permuted Observed Sample of Cell by sampling 100 times", ylab = "Pred observed values", xlab = "index")

hist(result$shdat0[,i], main = "Permuted Observed Sample", ylab = "Frequency", xlab = "Pred observed values")
qqnorm(result$shdat0[,i], main = "Permuted Observed Sample of cell 1 by sampling 100 times")
abline(a = 0, b = 1, lty = 2, col = "gray")

plot(result$shdat0[1,], main = "Observed Sample of all cells by 1st sampling", ylab = "Pred observed values", xlab = "index")




plot(result$rdat0[,i], main = "Rows Sample of one cell by sampling 100 times", ylab = "Pred observed values", xlab = "index")

hist(result$rdat0[,i], main = "Rows Sample", ylab = "Frequency", xlab = "Pred values")
qqnorm(result$rdat0[,i], main = "Rows Sample of cell 1 by sampling 100 times")

plot(result$rdat0[1,], main = "Rows Sample of all cells by 1st sampling", ylab = "Pred observed values", xlab = "index")



plot(result$cdat0[,i], main = "Columns Sample of one Cell by sampling 100 times", ylab = "Pred observed values", xlab = "index")

hist(result$cdat0[,i], main = "Columns Sample", ylab = "Frequency", xlab = "Pred values")
qqnorm(result$cdat0[,i], main = "Columns Sample of cell 1 by sampling 100 times")

plot(result$cdat0[1,], main = "Columns Sample: all cells by 1st sampling", ylab = "Pred observed values", xlab = "index")

i=1
plot(result$bdat0[,i], main = "Both rows and columns Sample \n of one cell by sampling 100 times", ylab = "Pred observed values", xlab = "index")

hist(result$bdat0[,i], main = "Rows and Columns Sample", ylab = "Frequency", xlab = "Pred values")

qqnorm(result$bdat0[,i], main = "Rows and Columns Sample of cell 1 by sampling 100 times")

plot(result$bdat0[1,], main = "Both rows and columns Sample: all cells by 1st sampling", ylab = "Pred observed values", xlab = "index")

# 
# 
# plot(result$shdat0)
# 
# 
# plot(c(result$shdat0[,1]))
# plot(result$rdat0[,1])
# plot(result$bdat0[, 1])
# plot(result$cdat0[, 1])
# 
# 
# plot(result$shdat0[1,])
# plot(result$rdat0[1,])
# plot(result$bdat0[1,])
# plot(result$cdat0[1,])
# 
# boxplot(result$shdat0[1,], result$rdat0[1,], result$cdat0[1,], result$bdat0[1,])
# boxplot(result$shdat0[,1], result$rdat0[,1], result$cdat0[,1], result$bdat0[,1])
# 
# 
# 
# for (i in seq(1,10, by=2)){
  boxplot(list(result$shdat0[, Asso_matrix[index] ==0], result$shdat0[,Asso_matrix[index] ==1]),
          main = paste("Permuted observed sample"),
          names = c("non-association", "association"),
          ylab = "Values",
          las = 1)
# }
# 
# 
  boxplot(list(result$shdat0[1,], result$shdat0[2,], result$shdat0[3,], result$shdat0[4,]),
          main = paste("Observed Values of all cells by a sampling"),
          names = c("1st per", "2nd per", "3rd per", "4th per"),
          ylab = "Observed Values",
          las = 2)
#   
#   #for (i in seq(1,10, by=2)){
  i=822
  i=389
  
  i=407
    boxplot(list(result$shdat0[, Asso_matrix[index] ==0], result$rdat0[, Asso_matrix[index] ==0], 
                 result$cdat0[, Asso_matrix[index] ==0], result$bdat0[, Asso_matrix[index] ==0]),
            main = paste("Permuted observed sample \n and Null distr for association = 0"),
            names = c("shdat0", "rdat0", "cdat0", "bdat0"),
            ylab = "Values",
            las = 1)
    
    
    boxplot(list(result$shdat0[, Asso_matrix[index] ==1], result$rdat0[, Asso_matrix[index] ==1], 
                 result$cdat0[, Asso_matrix[index] ==1], result$bdat0[, Asso_matrix[index] ==1]),
            main = paste("Permuted observed sample \n and Null distr for association = 1"),
            names = c("shdat0", "rdat0", "cdat0", "bdat0"),
            ylab = "Values",
            las = 1)
#   #}
#   
#     i=1
#     boxplot(list(result$rdat0[i,], result$cdat0[i,], result$bdat0[i,]),
#             main = paste("Null distr for all cells by 1st sampling"),
#             names = c("rdat0", "cdat0", "bdat0"),
#             ylab = "Values",
#             las = 1)
#   
# 
# plot(result$shdat0[,1], result$rdat0[,1])
# points(result$shdat0[true_values == 0], result$rdat0[true_values == 0], pch = 1, col = "blue")
# points(result$shdat0[true_values == 1], result$rdat0[true_values == 1], pch = 8, col = "red")
# 
# plot(result$shdat0)
# 
# 
# 
# 
# for (i in 1:4){
#   i=1
#   boxplot(list(result$shdat0[i,], result$rdat0[i,], result$cdat0[i,], result$bdat0[i,]),
#           main = paste("Pred Values:", paste0(i, "st"),"Sampling"),
#           names = c("shdat0", "rdat0", "cdat0", "bdat0"),
#           ylab = "Values",
#           las = 1)
# }
# 
# 
# 
# plot((result$shdat0[10,]), (result$rdat0[10,]))
# points(result$shdat0[true_values == 0], result$rdat0[true_values == 0], pch = 1, col = "blue")
# points(result$shdat0[true_values == 1], result$rdat0[true_values == 1], pch = 8, col = "red")
# 
# plot((result$shdat0[10,]), (result$cdat0[10,]))
# 
# plot((result$shdat0[10,]), (result$bdat0[10,]))
# 
# plot((result$rdat0[10,]), (result$bdat0[10,]))
# 



############### AUC - AUPR #######################

# #transpose to put the drugs in the column - If we do for leave-one-disease cross validation, this step is skipped
#pred_CV = prediction_matrix ## standard case

pred_CV = pred_CV_ComVal

inputObs_matrix=(Asso_matrix)
prediction_matrix=(pred_CV)

#
# #table(Observed=c(Asso_matrix)>0,Predicted=c(pred_CV)>0)
#
#
#sort inputObs_matrix by column using the decreasing order by column of prediction_matrix
res=sort_matrix(prediction_matrix,inputObs_matrix)
#res=sort_vector(prediction_matrix,inputObs_matrix)

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

AUC=trapz(fpr_list,tpr_list)#AUC

AUC


AUPR=trapz(recall_list,precision_list)#AUPR
AUPR


par(mfrow = c(1, 2))
plot(fpr_list,tpr_list, type="l", main=paste0("AUC=",round(AUC,3)), xlab="FPR (1-specificity)", ylab = "TPR (sensitivity)", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
abline(0,1,lty=2)

plot(recall_list,precision_list,type="l",main=paste0("AUPR=",round(AUPR,3)), xlab = "Recall", ylab="Precision", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
par(mfrow = c(1,1))

ks.omat = list()
ks.omat$fpr_list = fpr_list
ks.omat$tpr_list = tpr_list
ks.omat$recall_list = recall_list
ks.omat$precision_list = precision_list

## standard
  # NMF_standard_LRSSLdata = list()
  # NMF_standard_LRSSLdata$fpr_list = fpr_list
  # NMF_standard_LRSSLdata$tpr_list =          tpr_list
  # NMF_standard_LRSSLdata$recall_list    =       recall_list
  # NMF_standard_LRSSLdata$precision_list=           precision_list
  # save(NMF_standard_LRSSLdata, file = "Pred_NMF_standard_LRSSLdata.RData")
  # 
  
## Z-score
  # NMF_ttest_test2_Ydata = list()
  # NMF_ttest_test2_Ydata$fpr_list = fpr_list
  # NMF_ttest_test2_Ydata$tpr_list =          tpr_list
  # NMF_ttest_test2_Ydata$recall_list    =       recall_list
  # NMF_ttest_test2_Ydata$precision_list=           precision_list
  # save(NMF_ttest_test2_Ydata, file = "Pred_NMF_ttest_test2_Ydata.RData")
  
 ## Fisher
  # NMF_ttest_test3_LRSSLdata = list()
  # NMF_ttest_test3_LRSSLdata$fpr_list = fpr_list
  # NMF_ttest_test3_LRSSLdata$tpr_list =          tpr_list
  # NMF_ttest_test3_LRSSLdata$recall_list    =       recall_list
  # NMF_ttest_test3_LRSSLdata$precision_list=           precision_list
  # save(NMF_ttest_test3_LRSSLdata, file = "Pred_NMF_ttest_test3_LRSSLdata.RData")
  # 
 
##
  NMF_ttest_test4_hsdndata = list()
  NMF_ttest_test4_hsdndata$fpr_list = fpr_list
  NMF_ttest_test4_hsdndata$tpr_list =          tpr_list
  NMF_ttest_test4_hsdndata$recall_list    =       recall_list
  NMF_ttest_test4_hsdndata$precision_list=           precision_list
  #save(NMF_ttest_test4_hsdndata, file = "Pred_NMF_ttest_test4_hsdndata.RData")

  
  
  ### Load the fpt, tpr, recall, precision
  load("Pred_NMF_standard_oMat.RData")
  # load("Pred_NMF_standard_hsdn.RData")
  # load("Pred_NMF_standard_Cdata.RData")
  # load("Pred_NMF_standard_Fdata.RData")
  # load("Pred_NMF_standard_LRSSLdata.RData")
  # load("Pred_NMF_standard_Ydata.RData")
  
  standard = NMF_standard_oMat
  fpr_stan = standard$fpr_list
  tpr_stan = standard$tpr_list
  recall_stan = standard$recall_list
  precision_stan = standard$precision_list
  
  ### test 2 Z-score
  load("Pred_NMF_ttest_test2_oMatdata.RData")
  # load("Pred_NMF_ttest_test2_hsdndata.RData")
  # load("Pred_NMF_ttest_test2_Cdata.RData")
  # load("Pred_NMF_ttest_test2_Fdata.RData")
  # load("Pred_NMF_ttest_test2_LRSSLdata.RData")
  # load("Pred_NMF_ttest_test2_Ydata.RData")
  
  test2 = NMF_ttest_test2_oMatdata
  fpr_2 = test2$fpr_list
  tpr_2 = test2$tpr_list
  recall_2 = test2$recall_list
  precision_2 = test2$precision_list
  
  ## test 3, fisher
  load("Pred_NMF_ttest_test3_oMatdata.RData")
  # load("Pred_NMF_ttest_test3_hsdndata.RData")
  # load("Pred_NMF_ttest_test3_Cdata.RData")
  # load("Pred_NMF_ttest_test3_Fdata.RData")
  # load("Pred_NMF_ttest_test3_LRSSLdata.RData")
  # load("Pred_NMF_ttest_test3_Ydata.RData")
  
  test3 = NMF_ttest_test3_oMatdata
  fpr_3 = test3$fpr_list
  tpr_3 = test3$tpr_list
  recall_3 = test3$recall_list
  precision_3 = test3$precision_list
  
  
  ## test 4, 
  load("Pred_NMF_ttest_test4_oMatdata.RData")
  # load("Pred_NMF_ttest_test4_hsdndata.RData")
  
 
  test4 = NMF_ttest_test4_oMatdata
  fpr_4 = test4$fpr_list
  tpr_4 = test4$tpr_list
  recall_4 = test4$recall_list
  precision_4 = test4$precision_list
  
  
  
  ## bnnr
  load("omat_bnnr_test.RData")
  #load("hsdn_bnnr_test.RData")
  #load("Cdata_bnnr_test.RData")
  #load("Fdata_bnnr_test.RData")
  bnnr_test = omat_bnnr_test
  
  bnnr_fpr = as.matrix(bnnr_test$fpr)
  bnnr_tpr = as.matrix(bnnr_test$tpr)
  bnnr_recall = as.matrix(bnnr_test$recall)
  bnnr_precision = as.matrix(bnnr_test$precision)
  
  plot(fpr_stan, tpr_stan, type = "l", col = "red", lwd = 2, 
       xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC Curves")
  lines(fpr_2, tpr_2, col = "blue", lwd = 2)
  lines(fpr_3, tpr_3, col = "green", lwd = 2)
  #lines(fpr_4, tpr_4, col = "orange", lwd = 2)
  lines(bnnr_fpr, bnnr_tpr, col = "black", lwd = 2)
  lines(ks.omat$fpr_list, ks.omat$tpr_list, col = "orange", lwd=2)
  
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  plot.new()
  #legend("center", legend = c("Line 1", "Line 2"), col = c("blue", "red"), lty = 1, bty = 'n')
  legend("center", legend = c(paste("NMF standard: AUC =",round(trapz(fpr_stan,tpr_stan),2)),
                                   paste("NMF Per t.test with Z-score: AUC =",round(trapz(fpr_2,tpr_2),2)),
                                  paste("NMF Per t.test with Fisher's: AUC =",round(trapz(fpr_3,tpr_3),2)),
                                  #paste("NMF Per sum's: AUC =",round(trapz(fpr_4,tpr_4),2)),
                                  paste("NMF Per ks.test with Z-score: AUC =",round(trapz(ks.omat$fpr_list,ks.omat$tpr_list),2)),
                                  paste("BNNR: AUC =",round(trapz(bnnr_fpr,bnnr_tpr),2))),
                                        col = c("red", "blue", "green", "orange", "black"), lwd = 2)
  
  
  par(mfrow = c(1,1))
  
  plot(recall_stan, precision_stan, type = "l", col = "red", lwd = 2, 
       xlab = "False Positive Rate", ylim = c(0,0.4), ylab = "True Positive Rate", main = "Recall-Precision Curves")
  lines(recall_2, precision_2, col = "blue", lwd = 2)
  lines(recall_3, precision_3, col = "green", lwd = 2)
  #lines(recall_4, precision_4, col = "orange", lwd = 2)
  lines(ks.omat$recall_list, ks.omat$precision_list, col = "orange", lwd=2)
  lines(bnnr_recall, bnnr_precision, col = "black", lwd = 2)
  
  
  plot.new()
  #legend("center", legend = c("Line 1", "Line 2"), col = c("blue", "red"), lty = 1, bty = 'n')
  legend("center", legend = c(paste("NMF standard: AUPR =",round(trapz(recall_stan,precision_stan),2)),
                                   paste("NMF Per t.test with Z-score: AUPR =",round(trapz(recall_2,precision_2),2)),
                                   paste("NMF Per t.test with Fisher's: AUPR =",round(trapz(recall_3,precision_3),2)),
                                #paste("NMF Per sum: AUPR =",round(trapz(recall_4,precision_4),2)),
                              paste("NMF Per ks.test with Z-score: AUC =",round(trapz(ks.omat$recall_list,ks.omat$precision_list),2)),
                                paste("BNNR: AUPR =",round(trapz(bnnr_recall,bnnr_precision),2))),
         col = c("red", "blue", "green", "orange", "black"), lwd = 2)
  
  
  
  
  
##### qq plot

#qqnorm(prediction_matrix, main = "qq-Plot NMF oMat-MechDB dataset")
#qqnorm(prediction_matrix, main = "qq-Plot NMF hsdn-MechDB dataset")
#qqnorm(prediction_matrix, main = "qq-Plot NMF Cdataset")
#qqnorm(prediction_matrix, main = "qq-Plot NMF Fdataset")


# qqnorm(prediction_matrix, main = "qq-Plot NMF permuataion oMat-MechDB dataset")
# qqnorm(prediction_matrix, main = "qq-Plot NMF permuataion hsdn-MechDB dataset")
# qqnorm(prediction_matrix, main = "qq-Plot NMF permuataion Cdataset")
qqnorm(prediction_matrix, main = "qq-Plot NMF permuataion Fdataset")



### Brown's methods.


