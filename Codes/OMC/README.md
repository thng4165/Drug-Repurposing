# OMC


## Description: 
The original code was downloaded from https://github.com/BioinformaticsCSU/OMC/tree/master and includes the following functions:

*```Demo_OMC2.m``` demonstrates the experimental result on the gold standard dataset by OMC2 algorithm.

*```BNNR.m``` is the function of BNNR algorithm.

*```svt.m``` is the function of singular value thresholding operator.

*```KNN_drugS.m``` is the function of KNN Preprocessing based on drug similarity matrix.

*```KNN_diseaseS.m``` is the function of KNN Preprocessing based on disease similarity matrix.


We have added two additional functions:

*```sort_matrix.m```: Sorts or ranks the values in each column of a matrix.

*```TestBNNR_CV10folds.m```: Evaluates the OMC model using 10-fold cross-validation.

How to Run
Run "Demo_OMC2.m" to test the OMC algorithm on the sample dataset.

Run "TestOMC_CV10folds.m" to reproduce the 10-fold cross-validation results reported in our paper.
