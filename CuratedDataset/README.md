%the recent largest curated database from the Drug Mechanism Database project to collect drug-disease associations.
# Collection Databases
MechDB: drug - disease associations which are used to construct the drug-disease association matrix
DrugBank database: drugs and theirs structures which are used to calculate the drug similarity matrix
Orphanet database: disease - symptom associations of rare diseases which are used to construct the disease similarity matrix of oMat-MechDB data
Human Symptoms Disease Network (HSDN) database: disease - symptom associations of HSDN which are used to construct the disease similarity matrix of HSDN-MechDB

Please find these databases from https://drive.google.com/drive/folders/1ngBYPuepfTHPpbjvBpwHQ9G9-2NQJnb5?usp=sharing

# oMat-MechDB data
oMat-MechDB focuses on rare diseases with the disease symptoms collected from Orphanet database
Disease symptom dataset (version 2023) is collected from https://www.orphadata.com/alignments/

we collect the rare diseases and their symptom information from the Orphanet database (Version 2023) to create a symptom-disease matrix ($M_{sd}$). Next, we keep the diseases ($D$) that belong to both the drug-disease and symptom-disease sets. Finally, we retain drugs ($R$) that each has at least one association with the diseases in $D$. After filtering, the final association matrix comprise 89 diseases and 150 drugs, with 271 associations, 13079 non-associations, and the sparsity of 97.97\%.
The drug-drug similarity matrix $M_{rr}$ is calculated based on the SMILES structure using the Taminoto method.
Finally, we compute the disease-disease similarity matrix $M_{dd}$ based on the symptom-disease matrix ($M_{sd}$) using the Gaussian interaction profile (GIP) kernel approach \cite{van2011gaussian}.


# hsdn-MechDB data
HSDNMechDB uses the disease symptoms from the Human Symptoms Disease Network (HSDN) database [
For this dataset, we follow a similar approach as performed for the oMat-MechDB dataset, but construct a symptom-disease matrix ($M_{sd}$) using the Human Symptoms Disease Network (HSDN) database \cite{zhou2014human}, producing a large DR dataset.
The HSDN-MechDB dataset comprises 616 diseases and 1270 drugs with 3,710 associations and 778,619 non-associations, resulting in a high sparsity of 99.52\%, 
