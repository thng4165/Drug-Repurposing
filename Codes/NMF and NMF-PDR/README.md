# Non-negative Matrix Factorization (NMF) and Permutation-based NMF (NMF-PDR)
NMF is an ideal choice for drug repurposing models due to its ability to reveal hidden patterns and relationships within datasets, which is particularly important for uncovering unknown associations in the drug-disease association matrix. NMF aims to decompose a given non-negative matrix $V\in R^{n\times m}$ into the product of two non-negative matrices $W \in R^{n\times r}$ and $H \in R^{r\times m}$ such that $V \approx WH$, where rank $r$ is the dimension of drug feature and disease feature in the lower-rank spaces. This rank parameter determines the size of the new sub-matrices and is the primary parameter of NMF.

We introduce NMF-PDR, a proposed permutation approach for NMF in DR. NMF-PDR assumes that the drug-disease associations are not Evaluation of drug repurposing methods random, but are instead driven by biological and/or chemical relationships. Therefore, the prediction from NMF on the true drug-disease association matrix should be significantly higher than those from a random association matrix.
 
# Description

* ```drugRepurposing_functions.R```: Contains essential functions for the NMF and NMF-PDR models, including:

find_optimal_rank(): Determines the optimal rank for the NMF model using an SVD-based approach.

NMF_ELee(): Implements the NMF algorithm using Lee–Seung update rules.

sort_matrix(): Sorts or ranks the values in each column of a matrix.

makeFolds() and convert_element(): Utility functions for preparing 10-fold cross-validation.

* ```TestNMF_CV10folds.R```: Evaluates the NMF model using 10-fold cross-validation to assess performance

* ```TestNMF_PDR_CV10folds.R```: Evaluates the NMF-PDR model using 10-fold cross-validation to assess performance.

# How to Run

Run the appropriate script:

*```TestNMF_CV10folds.R``` for the NMF model

*```TestNMF_PDR_CV10folds.R``` for the NMF-PDR model

For each model, specify the input datasets, including:
* Drug–disease association matrix
* Disease similarity matrix
* Drug similarity matrix

For each model, the output including:
* A predicted score matrix
* Evaluation metrics: AUC and AUPR
