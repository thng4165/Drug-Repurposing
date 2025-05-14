# VDA_GMSBMF

# Drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization .
VDA_GMSBMF is a MF-based method which employs Gaussian kernel similarity and bilinear matrix factorization to explore potential virus-drug associations for SARS-CoV-2. Particularly, VDAGKSBMF applies Gaussian kernel similarity to the association matrix to enhance both virus and drug similarity, which improves the redictive capacity of bilinear matrix factorization [1]. This approach identifies new antiviral drugs by predicting unknown virus-drug associations and optimizing the model with the alternating-direction multiplier method (ADMM). Although originally designed for drug-virus associations, VDA-GKSBMF can also be applied to other drug disease association datasets, as used in this study.


## Description
The original code was downloaded from https://github.com/xiangju0208/VDA_GMSBMF/tree/main and includes the following functions:

"A_VDA_GMSBMF.m": The core implementation of the VDA_GMSBMF algorithm. 

"getKfoldCrossValidMatIndSet.m": 

"getPerfMetricROCcompute.m": 

"main.m": Performs cross-validation using the drug–virus dataset in VDdataset1.

We have added two additional functions:

"sort_matrix.m": Sorts or ranks the values in each column of a matrix.

"TestVDA_CV10folds.m": Evaluates the VDA_GMSBMF model using 10-fold cross-validation on drug-disease datasets.


## How to Run
Run "main.m" to test the VDA_GMSBMF algorithm using datasets in VDdataset1.

Run "TestVDA_CV10folds.m" to reproduce the 10-fold cross-validation results on the drug–disease dataset as reported in our paper.

## Codes 
A_VDA_GMSBMF.m: VDA-GKSBMF algorithm  <br>
% Input:  <br>
% matDV: drug-virus matrix <br> 
% Wdd: drug-drug matrix <br> 
% Wvv: virus-virus matrix <br> 
% others: parameters <br> 
% Ouput: <br>
% M_recovery: matrix for predicted drug-virus scores <br> 




## cite
[1] Yang, et al. Computational drug repositioning based on multi-similarities bilinear matrix factorization. Briefings in Bioinformatics 22.4 (2021): bbaa267. <br> 
Wang, et al. Drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization. <br>  

