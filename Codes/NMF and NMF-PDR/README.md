# Non-negative Matrix Factorization (NMF) and Permutation-based NMF (NMF-PDR)

## Description
*```drugRepurposing_functions.R```: Contains essential functions for the NMF and NMF-PDR models, including:

find_optimal_rank(): Determines the optimal rank for the NMF model using an SVD-based approach.

NMF_ELee(): Implements the NMF algorithm using Lee–Seung update rules.

sort_matrix(): Sorts or ranks the values in each column of a matrix.

makeFolds() and convert_element(): Utility functions for preparing 10-fold cross-validation.

*```TestNMF_CV10folds.R```: Evaluates the NMF model using 10-fold cross-validation to assess performance

*```TestNMF_PDR_CV10folds.R```: Evaluates the NMF-PDR model using 10-fold cross-validation to assess performance.

## How to Run
Run the appropriate script:
```TestNMF_CV10folds.R``` for the NMF model
```TestNMF_PDR_CV10folds.R``` for the NMF-PDR model

For each model, specify the input dataset, including:
The drug–disease association matrix
Disease similarity matrix
Drug similarity matrix

For each model, the output including
Each script generates:
A predicted score matrix
Evaluation metrics: AUC and AUPR
