setwd("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/oMatMechDB/DataPreparation/")

########################### drug-rare disease data##########################
## create the dataset including rare disease from oMat and drug from DrugMechDB such that
## diseases in DrugMechDB are inlcuded in oMat

load("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/chemInfo_subset_s.RData")
load("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/drug_disease_subset_s.RData")
load("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom/DRNMF/Data/InputDatasets/disease_symptom_oMat.RData")

rare_intersection = intersect(disease_symptom_oMat$mesh, drug_disease_subset_s$disease_mesh)
drug_rare_disease = drug_disease_subset_s[drug_disease_subset_s$disease_mesh %in% rare_intersection, ]

## Matrix indicating the interactions between drugs and diseases in subset drug_rare_disease (same disease in two datasets)
interact = table(drug_rare_disease$drugbank, drug_rare_disease$disease_mesh)
interact = as.matrix(interact)


## Matrix indicating the interaction between drugs and disease (-1/0/1 - unknow/no/yes interaction)
## Initial matrix  -1
in_drug_rare_disease_all = matrix(0, nrow = length(unique(drug_disease_subset_s$drugbank)), 
                                ncol = length(unique(disease_symptom_oMat$mesh))) ## 1319 * 1232
rownames(in_drug_rare_disease_all) = unique(drug_disease_subset_s$drugbank)
colnames(in_drug_rare_disease_all) = unique(disease_symptom_oMat$mesh)

## Combine interaction "intract" into the interaction matrix with all drugs and rare diseases.
in_drug_rare_disease_all[rownames(interact),colnames(interact)] = interact

heatmap(in_drug_rare_disease_all,
        col = colorRampPalette(c("white", "red"))(100),
        main = "Drug-Disease",
        xlab = "Disease",
        ylab = "Drug",
        margins = c(5,5),
        labRow = unique(drug_disease_subset_s$drugbank),
        labCol = unique(disease_symptom_oMat$mesh))

save(in_drug_rare_disease_all, file = "in_drug_rare_disease_all.RData")



