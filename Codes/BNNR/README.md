# BNNR
BNNR is a novel computational method, which utilizes Bounded Nuclear Norm Regularization algorithm to identify potential novel indications for known or new drugs. The code in this package implements Bounded Nuclear Norm Regularization (BNNR) for drug repositioning, which is implemented in Matlab2014a.


# Description

The original code was downloaded from https://github.com/BioinformaticsCSU/BNNR/tree/master and includes the following three functions:

"Demo.m": Demonstrates the experimental results on the Fdataset using the BNNR algorithm.

"BNNR.m": The core implementation of the BNNR algorithm.

"svt.m": Implements the Singular Value Thresholding (SVT) operator.

We have added two additional functions:

"sort_matrix.m": Sorts or ranks the values in each column of a matrix.

"TestBNNR_CV10folds.m":  Evaluates the BNNR model using 10-fold cross-validation.

# How to Run
Run "Demo.m" to test the BNNR algorithm on the sample dataset.

Run "TestBNNR_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.

