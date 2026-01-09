
# Collection Databases

* MechDB: Drug - disease associations which are used to construct the drug-disease association matrix. This is the recent largest curated database from the Drug Mechanism Database project to collect drug-disease associations.
* DrugBank database: drugs and their structures which are used to calculate the drug similarity matrix
* Orphanet database: disease - symptom associations of rare diseases which are used to construct the disease similarity matrix of oMat-MechDB data
* Human Symptoms Disease Network (HSDN) database: disease - symptom associations of HSDN which are used to construct the disease similarity matrix of HSDN-MechDB

Please find these databases from [https://drive.google.com/drive/folders/1ngBYPuepfTHPpbjvBpwHQ9G9-2NQJnb5?usp=sharing](https://drive.google.com/drive/folders/1Lo5Xg--8m2D-gbuzO9cDiaqxTfQa8Tlk?usp=sharing)

# oMat-MechDB data

oMat-MechDB focuses on rare diseases with the disease symptoms collected from Orphanet database. Disease-symptom dataset (version 2023) is collected from https://www.orphadata.com/alignments/. 

* Idea: We collect the rare diseases and their symptom information from the Orphanet database (Version 2023) to create a symptom-disease matrix ($M_{sd}$). Next, we keep the diseases ($D$) that belong to both the drug-disease and symptom-disease sets. Finally, we retain drugs ($R$) that each has at least one association with the diseases in $D$. After filtering, the final association matrix comprises 89 diseases and 150 drugs, with 271 associations, 13079 non-associations, and the sparsity of 97.97\%.
The drug-drug similarity matrix $M_{rr}$ is calculated based on the SMILES structure using the Taminoto method.
Finally, we compute the disease-disease similarity matrix $M_{dd}$ based on the symptom-disease matrix ($M_{sd}$) using the Gaussian interaction profile (GIP) kernel approach.

* Data preparation: The scripts ```Preparation_Database_oMat_MechDB.R``` and ```oMat_MechDB_processing.R``` prepare data, and construct the oMat-MechDB dataset, including the following components:
  - Drug–disease association matrix (```rare_disease_interact.RData```)  
  - Drug–drug similarity matrix (```rare_drug_sim.RData```)
  - Disease–disease similarity matrix (```rare_disease_sim.RData```)

* The folder ```Pre_datasets``` contains intermediate datasets generated during data processing. These datasets support the construction of the final oMat-MechDB data.

# HSDN-MechDB data

We follow a similar approach as performed for the oMat-MechDB dataset, but construct a symptom-disease matrix ($M_{sd}$) using the Human Symptoms Disease Network (HSDN) database [1]. The HSDN-MechDB dataset comprises 616 diseases and 1270 drugs with 3,710 associations and 778,619 non-associations, resulting in a high sparsity of 99.52\%, 

* Data preparation: The script ```hsdn_processing.R``` and ```Data_hsdn_MechDB.R``` prepares data, and constructs the hsdn-MechDB dataset, including the following components:
  - Drug–disease association matrix (```hsdn_MechDB_dd_association.RData```)  
  - Drug–drug similarity matrix (```hsdn_MechDB_drug_sim.RData```)
  - Disease–disease similarity matrix (```hsdn_MechDB_disease_sim_GIP.RData```)

* The folder ```Pre_datasets_hsdn``` contains intermediate datasets generated during data processing. These datasets support the construction of the final hsdn-MechDB data.

# References
[1] Zhou, XueZhong, Jörg Menche, Albert-László Barabási, and Amitabh Sharma. "Human symptoms–disease network." Nature communications 5, no. 1 (2014): 4212.

[2] Gonzalez-Cavazos, Adriana Carolina, Anna Tanska, Michael Mayers, Denise Carvalho-Silva, Brindha Sridharan, Patrick A. Rewers, Umasri Sankarlal, Lakshmanan Jagannathan, and Andrew I. Su. "DrugMechDB: a curated database of drug mechanisms." Scientific Data 10, no. 1 (2023): 632.
