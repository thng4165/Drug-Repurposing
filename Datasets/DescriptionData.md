# Description datasets
* Each model uses the following three input matrices:
    - Drug–disease Association Matrix
    - Disease–disease Similarity Matrix
    - Drug–drug Similarity Matrix
      
*The methods are implemented in both MATLAB and R. Accordingly, the input datasets are provided in both .mat (MATLAB) and .RData (R) formats.

## RDataFiles folder
### HDVD data
* ```HDVDdata_asso.RData```: drug-disease association matrix
* ```HDVDdata_disease_sim.RData```: disease-disease similarity matrix
* ```HDVDdata_drug_sim.RData```: drug-drug similarity matrix

### LAGCN data
* ```LAGCNdata_asso.RData```: drug-disease association matrix
* ```LAGCNdata_disease_sim.RData```: disease-disease similarity matrix
* ```LAGCNdata_drug_sim.RData```: drug-drug similarity matrix

### Fdata
* ```Fdata_asso.RData```: drug-disease association matrix
* ```Fdata_disease_sim.RData```: disease-disease similarity matrix
* ```Fdata_drug_sim.RData```: drug-drug similarity matrix

### Cdata
* ```Cdata_asso.RData```: drug-disease association matrix
* ```Cdata_disease_sim.RData```: disease-disease similarity matrix
* ```Cdata_drug_sim.RData```: drug-drug similarity matrix
  
### LRSSL
* ```LRSSLdata_asso.RData```: drug-disease association matrix
* ```LRSSLdata_disease_sim.RData```: disease-disease similarity matrix
* ```LRSSLdata_drug_sim.RData```: drug-drug similarity matrix

### Ydata
* ```Ydata_asso.RData```: drug-disease association matrix
* ```Ydata_disease_sim.RData```: disease-disease similarity matrix
* ```Ydata_drug_sim.RData```: drug-drug similarity matrix

### oMat-MechDB data
* ```rare_disease_interact.RData```: drug-disease association matrix
* ```rare_disease_sim.RData```: disease-disease similarity matrix
* ```rare_drug_sim.RData```: drug-drug similarity matrix
  
### HSDN-MechDB data
* ```hsdn_MechDB_dd_association.RData```: drug-disease association matrix
* ```hsdn_MechDB_disease_sim_GIP.RData```: disease-disease similarity matrix
* ```hsdn_MechDB_drug_sim.RData```: drug-drug similarity matrix

## MatlabDataFiles folder

### oMat-MechDB data
* ```rare_disease_interact.mat```: drug-disease association matrix
* ```rare_disease_sim.mat```: disease-disease similarity matrix
* ```rare_drug_sim.mat```: drug-drug similarity matrix

### HSDN-MechDB data
* ```hsdn_MechDB_dd_association_numeric.mat```: drug-disease association matrix
* ```hsdn_MechDB_disease_sim_GIP.mat```: disease-disease similarity matrix
* ```hsdn_MechDB_drug_sim.mat```: drug-drug similarity matrix

### Six public datasets
Six datasets HDVD (```HDVD.mat```), LAGCN (```LAGCN.mat```), Fdata (```Fdataset.mat```), Cdata (```Cdataset.mat```), LRSSL (```LRSSL.mat```), Ydata (```Ydataset.mat```) include:
* drug: drug-drug similarity matrix
* disease: disease-disease similarity matrix
* didr: drug-disease association matrix

## HGIMCdata
Three datasets containing multiple similarity measures for drugs and diseases: 

1. Cdataset_ms
2. Fdataset_ms
3. Ydataset_ms

Detail:

* Wrname: the DrugBank IDs of drugs;
* Wdname: the OMIM IDs of diseases;
* drug_ChemS: chemical structure similarity matrix;
* drug_AtcS: drug's ATC code similarity matrix;
* drug_SideS: side-effect similarity matrix;
* drug_DDIS: drug-drug interaction similarity matrix;
* drug_TargetS: drug's target profile similarity matrix;
* disease_PhS: disease phenotype similarity matrix;
* disease_DoS: disease ontology similarity matrix;
* didr: disease-drug association matrix.
