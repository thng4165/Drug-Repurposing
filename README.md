# Benchmarking heterogeneous network-based methods for drug repurposing

## Introduction
In this study, we systematically evaluate nine advanced heterogeneous network-based DR
methods across eight diverse datasets, including six publicly available datasets and two newly introduced
drug-disease datasets. The methods are based on (i) matrix factorization: NMF, NMF-PDR, NMF-DR
and VDA-GKSBMF, (ii) matrix completion: BNNR, OMC and HGIMC and (iii) recommendation system:
IBCF and LIBMF. We assess their performance using multiple evaluation metrics including the area under
the receiver operating characteristic (AUC) and area under the precision-recall curve (AUPR), analyze the
impact of data sparsity, and compare our findings with previous benchmarking studies.

## Description
This repository contains codes and datasets used in our study on benchmarking heterogeneous network-based methods for drug repurposing. It includes data preprocessing, model training, evaluation scripts, and result visualization.

### Codes/Methods
Each method is provided in a separate folder. A README file is included in each folder to describe the corresponding method and how to run it.


### Datasets

For each dataset, we provide data in two formats: MATLAB (.mat files) for methods implemented in MATLAB, and R (.RData files) for methods implemented in R.

The datasets include six public datasets (HDVD, LAGCN, Fdata, Cdata, Ydata, LRSSL) and two curated datasets (oMat-MechDB and HHSDN-MechDB).

Six public datasets were obtained from https://zenodo.org/records/8357512, where additional information is also available.

Details of two curated datasets are provided in the CuratedDataset folder.


## Evaluation Results



