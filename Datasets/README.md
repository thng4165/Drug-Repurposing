# LighGBM Prediction Models 
This repository contains code and data preparation scripts for training and evaluating lightGBM models to predict active compounds against the **WDR91** target, using various molecular fingerprints and similarity-based features.

## Data Preparation
### Training/Test Data for Each Fingerprint (FP) type
* ```WDR91_TrainData_drug_FP.R```:  Generates training dataset for each FP type.  
  The following files are produced: ```acEmb_WDR91_{FP}.RData```
* ```WDR91_Step1TestData_drug_FP.R```: Generates test dataset for each FP type.  
  The following files are produced: ```Test_WDR91_{FP}.RData```

### Calculation Similarity Matrix between Drugs and Ligands, Amino Acids, or Peptides
We calculate molecular similarities using various fingerprint (FP) types to construct features for the model, employing either the Tanimoto or cosine similarity methods (as specified by ```simType```).

* **Ligands**: Similarity matrix between 14 ligands and 9 FP types.  
  Rscript file: ```WDR91_similarity_drugs_ligands.R```  
  The following files are produced: ```{simType}_sim_ligands_dist_{FP}.RData```

* **Amino Acids**: Similarity matrix between 20 amino acids and 9 FP types.  
  Rscript file: ```WDR91_similarity_acids_fp.R```  
  The following files are produced: ```{simType}_sim_acids_{FP}.RData```

* **Peptides**: Similarity matrix between 605 unique peptides and 9 FP types.  
  Code in Rscript file: ```prepareData_set2.R```  
 
### Optimal Features for Models
Optimal feature was investigated for each case of feature type as below:
* FPs_only: ```FPtype="ECFP4"```  
* AAs_only: ```FPtype="AVALON"```  
* Ligands_only: ```FPtype="ECFP6"```  
* Peptides_only: ```FPtype = "AVALON"```
  
 The corresponding data preparation was implemented in ```prepareData_set2.R``` 

 ## Run LightGBM models
 ### Package
  ```lightgbm``` package and **R version 4.3.2** are used.
 ### Run Instructions
 Rscript file: ```Target2035_lightGBM_predictionModels.R```
 
 Following these steps:
* run ```source("../Target2035_functions.R")```: helper functions
* run ```source("prepareData_set2.R")```: data preparation
* Set the similarity method ```simType``` from ```simTypeList``` (we used simType = "tanimoto")
* Set the feature type ```fstype``` from ```featureSetting```
* Set training proportion ```p_use = 0.95```: divide data into a training set (95%) and a test set (5%)
* Set the number of cross-validation folds using ```foldNum``` (we used foldNum = 10)
* Cross-validation method: building block (BB) or standard random cross-validation by random
* Parameter settings of lightGBM
  We use ```learningrate = 0.1```, ```maxbin = 100```, ```numleaves = 75```, and ```forcecol_wise = TRUE```, otherwise, the default values are applied
* Run the prediction models with ```foldNum``` folds
