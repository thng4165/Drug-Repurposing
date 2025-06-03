# NMF-DR
NMF-DR improve the DR performance through optimizing the rank selection, initial values, prediction step of NMF. The method begins by building a
heterogeneous drug-disease association network, through the integration of drug and disease similarity networks. Different from the standard NMF, the NMFDR
methodology involves three key steps: (1) selecting a suitable factorization rank r using the minimum description length (MDL) criterion; (2) initializing
the factor matrices based on a SVD-based method; (3) predicting drug-disease association by combination of the NMF method with an accelerated hierarchical
alternating least squares (A-HALS) algorithm.

# Description

The original code was obtained from https://github.com/sshaghayeghs/NMF-DR/tree/main

We have added two additional functions:

* "sort_matrix.m": Sorts or ranks the values in each column of a matrix. It is used to calculate the AUC and AUPR for evaluation performance.

* "TestNMFDR_CV10folds.m": Evaluates the NMF-DR model using 10-fold cross-validation for performance assessment.

# How to Run
Run "TestNMFDR_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.
