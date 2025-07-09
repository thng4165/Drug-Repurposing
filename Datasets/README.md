# lighGBM models for 
This repository contains code and data preparation scripts for training and evaluating lightGBM models to predict active compounds against the **WDR91** target, using various molecular fingerprints and similarity-based features.
## Data preparation
### Training/Test Data for Each Fingerprint (FP)
* ```WDR91_TrainData_drug_FP.R```:  Generates training dataset for each fingerprint.
  The following files are produced: ```acEmb_WDR91_{FP}.RData```
* ```WDR91_Step1TestData_drug_FP.R```: Generates test dataset for each fingerprint.
  The following files are produced: ```Test_WDR91_{FP}.RData```

### Calculation Similarity Matrix between Drugs and Ligands, Amino Acids, or Peptides
We calculate molecular similarities using various fingerprint (FP) types to construct features for the model.
Ligands: Similarity between 14 ligands and 9 FP types.

Amino Acids: Similarity between 20 amino acids and 9 FP types.

Peptides: Similarity between 605 unique peptides and 9 FP types.

### Features for Models
Ligand (using ECFP6)
Amino Acid (using AVALON)
Peptide (using AVALON).
