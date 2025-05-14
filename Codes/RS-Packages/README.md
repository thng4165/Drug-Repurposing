# Recommender system packages
A recommendation system (RS) aims to predict user preferences and suggest relevant items based on historical data. It is widely used in e-commerce, streaming services, and social media. Drug repurposing and recommendation systems share a common mathematical foundation: predicting missing associations in a matrix.

In this study, we select two widely used recommendation system methods including item-based collaborative filtering (IBCF) and library for parallel matrix
factorization in shared memory systems (LIBMF) for application in drug repurposing.

## Description
*```drugRepurposing_functions.R```: Contains the necessary functions for testing.

*```Calculate_AUC_AUPR_func.R```: Function used to calculate AUC and AUPR.

*```TestRSpackages_CV10folds.R```: Evaluates the RS packages (ICBF and LIBMF) using 10-fold cross-validation. You can choose which package to run within the script.

## How to Run
Select the dataset and specify the package to run within the code.
Then, run "TestRSpackages_CV10folds.R" to reproduce the 10-fold cross-validation results reported in our paper.
