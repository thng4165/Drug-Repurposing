## in this R file, we claculate the similarity matrix of drug, extract the drug-disease dataset WELL-DONE


## load the dataset drug-disease
load("data/drug_disease_subset.RData")
load("data/chemInfo_subset.RData")



## intersection of disease between disease-symptom-oMat, drug-disease-subset (small dataset)
## We did again in RareDiseaseDrug.R
intersect_disease = intersect(unique(disease_symptom_oMat$mesh), unique(drug_disease_subset$disease_mesh))
length(intersect_disease)
## datasets with the same diseases in case we would like to see the interaction
extract_disease_symptom = disease_symptom_oMat[disease_symptom_oMat$mesh %in% intersect_disease, ]
extract_drug_disease = drug_disease_subset[drug_disease_subset$disease_mesh %in% intersect_disease, ]
extract_chemInfo = chemInfo_subset[chemInfo_subset$DrugBankID %in% extract_drug_disease$drugbank, ]


####################################### Drug-Drug similarity Tanimoto for all drug in drug-disese after extracting drugbank, removing non-smiles, .....

#install.packages("BiocManager")
#BiocManager::install("ChemmineOB")
#install.packages("ChemmineR")
##install.packages("Cheminformatics")

library(ChemmineR)
#library(Cheminformatics)


mysmiles = chemInfo_subset$SMILES
names(mysmiles) = chemInfo_subset$DrugBankID

sdfset = smiles2sdf(mysmiles) ### ==> 13 invalid SDFs detected ==> remove
valid = validSDF(sdfset)  ## check the valide entries
sdfset = sdfset[valid] ## only keep the valide entries and remove the invalide entries
apset = sdf2ap(sdfset) 
fpset = desc2fp(apset)
drug_sim = sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE)) ## drug-drug similarity
dim(drug_sim) ## 1319 x 1319 (remove 13 invalid SDFs strucutures)
drug_sim[1:5,1:5]

heatmap(drug_sim,
        col = colorRampPalette(c("white", "red"))(100),
        main = "Drug - Drug",
        xlab = "Drugs",
        ylab = "Drug",
        margins = c(5, 5))

save(drug_sim, file = "data/drug_sim.RData") ## WELL DONE

############################## Extract sub-data again ###########################
## Since 13 SMILES of 13 drugs have invalid SDF structures, and we remove them. it means that the chemInfo_subset
## is reduced 13 drugs. We extract again the drug_disease_subset then two these datasets have the same list of drugs

invalid_indices = which(!valid) ## indices of invalid term

## dataset after removing drug with invalid SMILES --> SDF
chemInfo_subset_s = chemInfo_subset[-invalid_indices, ] ## drugbank data after remove 
removed_DrugBank = chemInfo_subset$DrugBankID[invalid_indices] ## indicates of removed elements
drug_disease_subset_s = drug_disease_subset[!drug_disease_subset$drugbank %in% removed_DrugBank, ]
dim(drug_disease_subset_s) ## 4124 x 11 

save(chemInfo_subset_s, file = "data/chemInfo_subset_s.RData") # after remove invalid validSDF
save(drug_disease_subset_s, file="data/drug_disease_subset_s.RData") # after remove invalid validSDF
## WELL-DONE

################################ Interaction between disease and drug #############################
# ## Matrix indicating the interactions between drugs and diseases from DB data, we does not use....
# interaction_DrugDisease_subset <- table( drug_disease_subset_s$disease_mesh, drug_disease_subset_s$drugbank)
# interaction_DrugDisease_subset <- as.matrix(interaction_DrugDisease_subset)
# dim(interaction_DrugDisease) ## 672 * 1319 (unique(disease) * unique(drug))
# dim(drug_disease_subset)
# heatmap(interaction_DrugDisease_subset,
#         col = colorRampPalette(c("white", "red"))(100),
#         main = "Drug-Disease",
#         xlab = "Drug",
#         ylab = "Disease",
#         margins = c(5, 5),
#         labRow = unique(drug_disease_subset$disease_mesh),
#         labCol = unique(drug_disease_subset$drugbank))
# 
# save(interaction_DrugDisease_subset, file = "interaction_DrugDisease_subset.RData")
# 


