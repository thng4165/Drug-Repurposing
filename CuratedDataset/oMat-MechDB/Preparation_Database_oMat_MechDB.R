### collect and build a network between drug-disease-symptom for both oMat-MechDB and hsdn-MechDB data
rm(list=ls())
setwd("C:/Trang/KIProjects/RareDiseases/oMatMechDB1/oMatMechDB/DataPreparation")


##---------------- Orphanet database: oMat ---------------------------------
## All data can be found here: https://www.orphadata.com/orphanet-scientific-knowledge-files/
## 16 Aug 2023: Download disease-symptom data en_product4.xml from https://www.orphadata.com/phenotypes/
## OrphaNet data contains databases for disease-symptom network, focusing mostly on rare diseases 
## Creating "orphaNet.RData" and "orphaNet_Disease_mapping.RData" including disease-symptom and their information for rare diseases data
 
##---------------------Processing with "en_product4.xml"
library("xml2")
dbdb=read_xml("Data/Orphanet/en_product4.xml")
kids <- xml_children(dbdb)
length(kids)
d=2
x=as_list(kids[[d]])
head(x)

orphaNet_db=x
oMat=NULL
for (i in 1:length(x)){
  if (i %% 100 == 0) cat (" ",i)
  #i=1
  xi=x[[i]]
  #names(xi)
  y=xi$Disorder
  #names(y)
  diseaseName=y$Name[[1]]
  OrphaCode=as.integer(y$OrphaCode[[1]])
  phenNum=as.integer(attributes(y$HPODisorderAssociationList)$count)
  if (phenNum > 0){
    for (k in 1:phenNum){
      #k=1
      y1=y$HPODisorderAssociationList[[k]]
      attributes(y1)
      names(y1$HPO)
      phen=data.frame(OrphaCode=OrphaCode,Disease=diseaseName,HPOId=y1$HPO$HPOId[[1]],HPOTerm=y1$HPO$HPOTerm[[1]],Frequency=y1$HPOFrequency$Name[[1]])
      oMat=rbind(oMat,phen)
    }
  }
}

# save(orphaNet_db, oMat,file="Pres_datasets/orphaNet.RData")


##---------------- Processing with "en_product1.xml"
## orphaNet data contains disease mapping with 7 other databases including:
## ICD-10, ICD-11, OMIM, UMLS, MeSH, MedDRA and GARD
library("xml2")
dbdb=read_xml("Data/Orphanet/en_product1.xml")
kids <- xml_children(dbdb)
length(kids)
d=2
x=as_list(kids[[d]])
orphaNet_diseaseMap=x

orphaNet_OrphaCode=sapply(orphaNet_diseaseMap, function(x) unlist(x$OrphaCode))
names(orphaNet_OrphaCode)=NULL
which(orphaNet_OrphaCode==586)


# library("rjson")
# orphaNet_diseaseMap <- fromJSON(file="Data/Orphanet/en_product1.json")
# orphaNet_diseaseMap=orphaNet_diseaseMap$JDBOR[[1]]
# names(orphaNet_diseaseMap)

## now we get all synonyms of a diease
orphaNet_synonyms=lapply(orphaNet_diseaseMap, function(x){
  x1=trimws(c(unlist(x$Name),unlist(x$SynonymList)))
  names(x1)=NULL
  x1=x1[x1!=""]
  return(x1)
})
head(orphaNet_synonyms)
length(orphaNet_synonyms)
names(orphaNet_synonyms)=orphaNet_OrphaCode


## now we get MeSH and other IDs of a diease
orphaNet_allIDs=lapply(orphaNet_diseaseMap, function(x){
  dbMap=NULL
  OrphaCode=unlist(x$OrphaCode)
  #  cat("\n",OrphaCode)
  x1=x$ExternalReferenceList
  if (length(x1)>0){
    isExist=FALSE
    for (i in 1:length(x1)){
      x2=x1[[i]]
      if (length(grep("Source",names(x2)))>0){
        s=x2$Source[[1]]
        id=x2$Reference[[1]]
        dbMap=rbind(dbMap,c(s,id))
        isExist=TRUE
      }
    }
    if (isExist){
      colnames(dbMap)=c("DB","ID")
      dbMap=as.data.frame(dbMap)
      rownames(dbMap)=NULL
    }
  }
  return(dbMap)
})
head(orphaNet_allIDs)
length(orphaNet_allIDs)
names(orphaNet_allIDs)=orphaNet_OrphaCode

# save(orphaNet_allIDs,orphaNet_synonyms, file="orphaNet_Disease_mapping.RData")


##------------------------Disease-symptom
library("ontologyIndex")
library("ontologySimilarity")
library("ontologyPlot")
file="Data/Orphanet/ORDO_en_4.3.obo"
orphanet_ORDO <- get_ontology(file)
names(orphanet_ORDO)
length(orphanet_ORDO$name)
p=which(orphanet_ORDO$id=="Orphanet:586")
p
orphanet_ORDO$name[p]

file="Data/Orphanet/HOOM_en_2.0.obo"
orphanet_HOOM <- get_ontology(file)
names(orphanet_HOOM)


# library("rjson")
# hpDat_raw <- fromJSON(file="Data/hp.json")
# hpo_onto=hpDat_raw$graphs[[1]]
# names(hpo_onto)


file="Data/hp.obo"
hpo_onto = get_ontology(file)
names(hpo_onto)

library(data.table)
hpo_anno <- fread(file="Data/phenotype.hpoa", skip = 4)
head(hpo_anno)
dim(hpo_anno)

# save(hpo_onto, hpo_anno,file="HPO.RData")

##---------------------- DrugMechBD database--------------------------
## DrugMechDB curated the drug-disease associations (https://www.nature.com/articles/s41597-023-02534-z)

##--------------- DrugMechDB data
library("rjson")
DrugMechDB <- fromJSON(file="Data//DrugMechDB-2.0.1//indication_paths.json")
length(DrugMechDB)
x=DrugMechDB[[999]]
names(x)
x$nodes


##-------------- DrugBank data
library("xml2")
path="Data/DrugBank/full_database_04Jan2023.xml"

dbdb=read_xml(path)
DrugBank=xml_contents(dbdb)

getIndication <-function(d){
  if (d %% 100 == 0) cat (" ",d)
  x=DrugBank[[d]]
  x1=xml_contents(x)
  
  indicationOut=nameOut=descriptionOut=NA
  
  p=grep("<indication",x1[1:(length(x1)-30)])
  x2=xml_contents(x1[p])
  x3=as_list(x2)
  if (length(x3)>0) indicationOut=x3[[1]][1]
  
  p=grep("<name",x1[1:(length(x1)-30)])
  x2=xml_contents(x1[p])
  x3=as_list(x2)
  if (length(x3)>0) nameOut=x3[[1]][1]
  
  p=grep("<description",x1[1:(length(x1)-30)])
  x2=xml_contents(x1[p])
  x3=as_list(x2)
  if (length(x3)>0) descriptionOut=x3[[1]][1]
  #  nameOut=NA
  #  descriptionOut=NA
  
  return(c(nameOut,descriptionOut, indicationOut))
}

indicationDbList=lapply(1:length(DrugBank),getIndication)
indicationDb=do.call(rbind,indicationDbList)
dim(indicationDb)
colnames(indicationDb)=c("DrugName","Description","Indication")

#save(DrugBank, indicationDb, file="DrugBank_drug_indication.RData")

##-------------- now extract disease names from the indication
x1=toupper(indicationDb[1,]$Indication)
x1=toupper(indicationDb[2,]$Indication)
x1
y1=toupper(unique(oMat$HPOTerm))
y2=toupper(unique(oMat$Disease))

z1=sapply(y1,function(y) grep(y,x1))
z2=lengths(z1)
p=which(z2>0)
length(p)
y1[p]

z1=sapply(y2,function(y) grep(y,x1))
z2=lengths(z1)
p=which(z2>0)
length(p)
y2[p]

HPOTerm=HPODisease=list()
for (i in 1:nrow(indicationDb)){ #15235 drugs
  if (i %% 100 == 0) cat (" ",i)
  x1=toupper(indicationDb[i,]$Indication)
  
  z1=sapply(y1,function(y) grep(y,x1))
  z2=lengths(z1)
  p=which(z2>0)
  #length(p)
  HPOTerm[[i]]=y1[p]
  
  z1=sapply(y2,function(y) grep(y,x1))
  z2=lengths(z1)
  p=which(z2>0)
  #length(p)
  HPODisease[[i]]=y2[p]
  
}

DrugBank_HPODisease=HPODisease
DrugBank_HPOTerm=HPOTerm
save(DrugBank, indicationDb, DrugBank_HPODisease, DrugBank_HPOTerm, file="DrugBank_drug_indication_mapping.RData")

# 
# ##---------------------------MESH databases--------------------------------
# ## Download mesh data version 2023: https://www.nlm.nih.gov/databases/download/mesh.html
# path="Data/desc2023.gz"
# dbdb=read_xml(path)
# 
# Mesh_xml=xml_contents(dbdb)
# # length(Mesh_xml)
# # # Mesh_xml=xml_contents(dbdb) #this is faster
# # # length(Mesh_xml)
# # # d=1
# # # x=Mesh_xml[[d]]
# # # x1=xml_contents(x)
# # # p=6
# # # x2=xml_contents(x1[p])
# # # x3=xml_contents(x2[1])
# # # x3
# 
# Mesh <- xml_children(dbdb) #this is slower but more structured
# # length(Mesh) # 30454 cases
# # d=1000
# # x=as_list(Mesh[[d]])
# # names(x)
# # x$DescriptorUI #id
# # unlist(x$DescriptorName) #name of drug/disease
# # #x$ConceptList
# 
# Mesh_id=sapply(Mesh_xml,function(x){
#   x1=xml_contents(x)
#   x2=xml_contents(x1[1])#DescriptorUI
#   x3=as.character(x2)
#   return(x3)
# })
# 
# Mesh_des=sapply(Mesh_xml,function(x){
#   x1=xml_contents(x)
#   x2=xml_contents(x1[2])#DescriptorName
#   x3=as.character(xml_contents(x2))
#   return(x3)
# })
# 
# Mesh_entryTerms=sapply(Mesh_xml,function(x){
#   x1=xml_contents(x)
#   x2=xml_contents(x1[2])#DescriptorName
#   x3=as.character(xml_contents(x2))
#   return(x3)
# })
# 
# 
# # save(Mesh_entryTerms, Mesh_des, Mesh_id, Mesh, Mesh_xml, file="Mesh.RData")
# 
# 
# ##-----------------------Human symptoms-disease network (hsdn) database----------------------------
# ## Human symptoms-disease network https://www.nature.com/articles/ncomms5212
# ## This database focuses on cancer disease, using Mesh id
# library(data.table)
# hsdn_net=fread(file="network_Human_symptoms–disease.txt")
# hsdn_disease=fread("diseases_Human_symptoms–disease.txt")
# hsdn_symptom=fread("symptoms_Human_symptoms–disease.txt")
# 
# 
# head(hsdn_net)
# head(hsdn_disease)
# head(hsdn_symptom)
# 
# # save(hsdn_net,hsdn_disease, hsdn_symptom, file="hsdn.RData")


