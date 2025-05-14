# NMF-DR

## Motivation:
Drug repurposing is a potential alternative to the traditional drug discovery process. Drug repurposing can be formulated as a recommender system that recommends novel indications for available drugs based on known drug-disease associations. This paper presents a method based on non-negative matrix factorization (NMF-DR) to predict the drug-related candidate disease indications. This work proposes a recommender system-based method for drug repurposing (NMF-DR) to predict novel drug indications by integrating drug and diseases related data sources.For this purpose, this framework first integrates two types of disease similarities, the associations between drugs and diseases, and the various similarities between drugs from different views to make a heterogeneous drug-disease interaction network. Then, an improved non-negative matrix factorization-based method is proposed to complete the drug-disease adjacency matrix with predicted scores for unknown drug-disease pairs.

## Description

The original code was downloaded from https://github.com/sshaghayeghs/NMF-DR/tree/main

We have added two additional functions:

"sort_matrix.m": Sorts or ranks the values in each column of a matrix. It is used to calculate the AUC and AUPR for evaluation performance.

"TestNMFDR_CV10folds.m": Evaluates the NMF-DR model using 10-fold cross-validation for performance assessment.

How to Run
Run "TestNMFDR_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.
