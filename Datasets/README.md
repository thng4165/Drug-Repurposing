# lighGBM models for 
This repository contains code and data preparation scripts for training and evaluating lightGBM models to predict active compounds against the **WDR91** target, using various molecular fingerprints and similarity-based features.
## Data preparation
### Training/Test Data for Each Fingerprint (FP) type
* ```WDR91_TrainData_drug_FP.R```:  Generates training dataset for each FP type.  
  The following files are produced: ```acEmb_WDR91_{FP}.RData```
* ```WDR91_Step1TestData_drug_FP.R```: Generates test dataset for each FP type.  
  The following files are produced: ```Test_WDR91_{FP}.RData```

### Calculation Similarity Matrix between Drugs and Ligands, Amino Acids, or Peptides
We calculate molecular similarities using various fingerprint (FP) types to construct features for the model.

* Ligands: Similarity matrix between 14 ligands and 9 FP types.  
  Rscript file: ```WDR91_similarity_drugs_ligands.R```  
  The following files are produced: ```{simType}_sim_ligands_dist_{FP}.RData```

* Amino Acids: Similarity matrix between 20 amino acids and 9 FP types.  
  Rscript file: ```WDR91_similarity_acids_fp.R```  
  The following files are produced: ```{simType}_sim_acids_{FP}.RData```

* Peptides: Similarity matrix between 605 unique peptides and 9 FP types.  
  Code in Rscript file: ```prepareData_set2.R```  
  The following files are produced:

### Optimal Features for Models
Optimal feature is investigated for each case as below:
* FPs_only: ```FPtype="ECFP4"```  
* AAs_only: ```FPtype="AVALON"```  
* Ligands_only: ```FPtype="ECFP6"```  
* Peptides_only: ```FPtype = "AVALON"```
  
 The corresponding data preparation is implemented in ```prepareData_set2.R``` 
