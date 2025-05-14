# VDA_GMSBMF

# Drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization .
VDA_GMSBMF is a MF-based method which employs Gaussian kernel similarity and bilinear matrix factorization to explore potential virus-drug associations for SARS-CoV-2. Particularly, VDAGKSBMF applies Gaussian kernel similarity to the association matrix to enhance both virus and drug similarity, which improves the redictive capacity of bilinear matrix factorization [1]. This approach identifies new antiviral drugs by predicting unknown virus-drug associations and optimizing the model with the alternating-direction multiplier method (ADMM). Although originally designed for drug-virus associations, VDA-GKSBMF can also be applied to other drug disease association datasets, as used in this study.


## Description
The original code was downloaded from https://github.com/BioinformaticsCSU/BNNR/tree/master and includes the following three functions:

"Demo.m": Demonstrates the experimental results on the Fdataset using the BNNR algorithm.

"BNNR.m": The core implementation of the BNNR algorithm.

"svt.m": Implements the Singular Value Thresholding (SVT) operator.

We have added two additional functions:

"sort_matrix.m": Sorts or ranks the values in each column of a matrix.

"TestBNNR_CV10folds.m": Evaluates the BNNR model using 10-fold cross-validation.

## How to Run
Run "Demo.m" to test the BNNR algorithm on the sample dataset.

Run "TestBNNR_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.

## Codes 
A_VDA_GMSBMF.m: VDA-GKSBMF algorithm  <br>
% Input:  <br>
% matDV: drug-virus matrix <br> 
% Wdd: drug-drug matrix <br> 
% Wvv: virus-virus matrix <br> 
% others: parameters <br> 
% Ouput: <br>
% M_recovery: matrix for predicted drug-virus scores <br> 

main.m: cross-validation code.  <br>


## Dataset
Data is located in the directory: ./VDdataset1 <br>    
matDrugVirus.txt <-- drug-virus matrix   <br> 
matDrugDrug.txt <-- drug-drug matrix <br> 
matVirusVirus.txt <-- virus-virus matrix   <br> 


## Results 
The results will be automatically saved into the directory: Results.   <br>

## cite
[1] Yang, et al. Computational drug repositioning based on multi-similarities bilinear matrix factorization. Briefings in Bioinformatics 22.4 (2021): bbaa267. <br> 
Wang, et al. Drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization. <br>  

