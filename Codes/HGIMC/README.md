# HGIMC

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

How to Run
Run "Demo_HGIMC.m" to test the HGIMC algorithm on the sample dataset.

Run "TestHGIMC_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.
