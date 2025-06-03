# HGIMC
HGIMC also utilizes BNNR, but combines with the guilt-by association principle of HGBI. It refines drug and disease similarity
matrices using Gaussian radial basis functions (GRB) before applying bounded matrix completion (BMC) with the optimization equation of
BNNR model to impute high confidence drug-disease associations. This step enriches the edges connecting drug and disease networks. Finally, it integrates the
updated drug and disease similarity matrices with the updated drug-disease association matrix to predict the unknown associations.

Originally, HGIMC was designed to integrate multiple similarity measures for drugs and diseases. However, since this study focuses on methods using a
single similarity measure, we apply HGIMC in a single similarity mode. Our analysis also shows minimal performance differences between using single and
multiple similarity measures.

# Description

The original code was downloaded from https://github.com/BioinformaticsCSU/HGIMC and includes the following functions:

* ```fBMC.m```: this function can implement the bounded matrix completion algorithm;
* ```fGRB.m```: this function can implement the Gaussian radial basis function;
* ```fHGI.m```: this function can implement the heterogeneous graph inference algorithm;
* ```fNorm.m```: this function can normalize the similarity matrix;
* ```svt.m```: this function can implement singular value thresholding operator.
* ```Demo_HGIMC.m```: Demonstrates the experimental results on the Fdataset_ms using the HGIMC algorithm.

We have added two additional functions:

* ```sort_matrix.m```: Sorts or ranks the values in each column of a matrix.

* ```TestBNNR_CV10folds.m```: Evaluates the HGIMC model for using 10-fold cross-validation.

# How to Run
* Run "Demo_HGIMC.m" to test the HGIMC algorithm on the sample dataset.

* Run "TestHGIMC_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.
