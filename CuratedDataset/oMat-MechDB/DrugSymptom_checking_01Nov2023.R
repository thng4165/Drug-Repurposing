### collect and build a network between drug-disease-symptom
rm(list=ls())
#setwd("C:/Nghiavtr/Workplaces/Projects/DrugRepurposing/DrugSymptom")
# get disease-symptom of rare data
load("orphaNet.RData")
head(oMat)

#get information of drugs their indications from DrugBank
load("DrugBank_drug_indication_mapping.RData")
head(indicationDb)

### other data of symptom-disease    
# Human symptoms-disease network https://www.nature.com/articles/ncomms5212
#this database focuses on cancer disease, using Mesh id
load("hsdn.RData")
head(hsdn_net)
head(hsdn_disease)
head(hsdn_symptom)


# other Disease-drug association network: #https://snap.stanford.edu/biodata/datasets/10004/10004-DCh-Miner.html


############################### 
### see codes for parsing the databases below, but so messy, sorry

######### all data can be found here: https://www.orphadata.com/orphanet-scientific-knowledge-files/
##### 16 Aug 2023: Download disease-symptom data en_product4.xml from https://www.orphadata.com/phenotypes/
### orphaNet data contains databases for disease-symptom network, focusing mostly on rare diseases 
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


save(orphaNet_db, oMat,file="orphaNet.RData")


### orphaNet data contains disease mapping with 7 other databases including ICD-10, ICD-11, OMIM, UMLS, MeSH, MedDRA and GARD
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

length(x)
x1=x[[82]]
names(x1)
x1$Name
x1$SynonymList
x2=x1$ExternalReferenceList
names(x2)
length(x2)
x3=x2[[5]] #corresponding disease from MeSH
names(x3)
x3$Source
x3$Reference
x3
unlist(sapply(x2, function(x) x$Source))

p=which(orphanet_ORDO$id=="Orphanet:586") #"Cystic fibrosis"
orphanet_ORDO$name[p]

# library("rjson")
# orphaNet_diseaseMap <- fromJSON(file="Data/Orphanet/en_product1.json")
# orphaNet_diseaseMap=orphaNet_diseaseMap$JDBOR[[1]]
# names(orphaNet_diseaseMap)

### now we get all synonyms of a diease
orphaNet_synonyms=lapply(orphaNet_diseaseMap, function(x){
  x1=trimws(c(unlist(x$Name),unlist(x$SynonymList)))
  names(x1)=NULL
  x1=x1[x1!=""]
  return(x1)
  })
head(orphaNet_synonyms)
length(orphaNet_synonyms)
names(orphaNet_synonyms)=orphaNet_OrphaCode



### now we get MeSH and other IDs of a diease
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


save(orphaNet_allIDs,orphaNet_synonyms, file="orphaNet_Disease_mapping.RData")

# #################################### use ROBOT in UPPMAX to convert owl to obo
# #### download data from orphanet: https://www.orphadata.com/ordo/ , phenotypes https://www.orphadata.com/phenotypes/ and hoom https://www.orphadata.com/hoom/
# 
# # we need to convert owl to obo files using ROBOT to be able to read in R
# cd /proj/snic2020-6-4/Nghia/Working/Xcelerate_RARE/ROBOT
# wget https://github.com/ontodev/robot/releases/download/v1.9.4/robot.jar
# curl https://raw.githubusercontent.com/ontodev/robot/master/bin/robot > robot
# chmod u+x robot
# 
# cd /proj/snic2020-6-4/Nghia/Working/Xcelerate_RARE/Orphanet
# wget https://www.orphadata.com/data/xml/en_product4.xml
# wget https://www.orphadata.com/data/ontologies/ordo/last_version/ORDO_en_4.3.owl
# wget https://www.orphadata.com/data/ontologies/hoom/hoom_orphanet_2.0.zip
# unzip hoom_orphanet_2.0.zip
# 
# #do conversion
# module load java/OracleJDK_11.0.9 #use more recent Java for the owl files
# #use --check false option to ignore the errors of obo structure error
# ../ROBOT/robot convert --input ORDO_en_4.3.owl --output ORDO_en_4.3.obo --check false
# ../ROBOT/robot convert --input HOOM_en_2.0.owl --output HOOM_en_2.0.obo --check false 

####################################
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



##### Look up parent HPO phenotypes to get BROAD categories: hp.json and phenotype.hpoa (July 2023)
#https://hpo.jax.org/app/data/ontology
#https://hpo.jax.org/app/data/annotations

# library("rjson")
# hpDat_raw <- fromJSON(file="Data/hp.json")
# hpo_onto=hpDat_raw$graphs[[1]]
# names(hpo_onto)

#install.packages("ontologyIndex") #instruction here: https://cran.r-project.org/web/packages/ontologyIndex/vignettes/reading-ontologies.html
#install.packages(c("ontologyIndex","ontologySimilarity","ontologyPlot"))

library("ontologyIndex")
library("ontologySimilarity")
library("ontologyPlot")

file="Data/hp.obo"
hpo_onto <- get_ontology(file)
names(hpo_onto)

library(data.table)
hpo_anno <- fread(file="Data/phenotype.hpoa", skip = 4)
head(hpo_anno)
dim(hpo_anno)


save(hpo_onto, hpo_anno,file="HPO.RData")

# ### do some plots with hpo data
 library("ontologyIndex")
 library("ontologySimilarity")
 library("ontologyPlot")
setwd("D:/TRANG/KI Projects/DrugSymptom/DrugSymptom")
load("Data/HPO.Rdata")
# #compute semetic similarity
 term_sets <- replicate(simplify=FALSE, n=5, expr=minimal_set(hpo_onto, sample(hpo_onto$id, size=8)))
 sim_mat <- get_sim_grid(ontology=hpo_onto, term_sets=term_sets)
# 
# #now explore some plots
# term_vec=unlist(term_sets)
# onto_plot(hpo_onto,terms=term_vec[1:2])
# onto_plot(hpo_onto,terms=get_ancestors(hpo_onto,term_vec[1:2]))
# 
# onto_plot(hpo_onto, terms=remove_links(hpo_onto, get_ancestors(hpo_onto, term_vec[1:2])))
# onto_plot(hpo_onto, terms=remove_links(hpo_onto, get_ancestors(hpo_onto, term_vec[1:5])))
# 
# #consider phenotypes of 3 diseases
# hpo_phenotypes <- list(
#   A=c("HP:0001382","HP:0004272","HP:0007917","HP:0004912","HP:0001596"),
#   B=c("HP:0001382","HP:0004272","HP:0002165","HP:0004800","HP:0004912"),
#   C=c("HP:0004800","HP:0001382","HP:0004912","HP:0007917","HP:0008743")
# )
# 
# onto_plot(hpo_onto, term_sets=hpo_phenotypes, label=label_by_term_set)
# onto_plot(
#   hpo_onto,
#   frequencies=get_term_frequencies(hpo_onto, hpo_phenotypes),
#   term_sets=hpo_phenotypes,
#   label=label_by_term_set,
#   edge_attributes=list(color="red", lty="dashed"))
# 
# 
# 

#### drug-disease network from DrugMechDB - https://www.nature.com/articles/s41597-023-02534-z
library("rjson")
DrugMechDB <- fromJSON(file="Data//DrugMechDB-2.0.1//indication_paths.json")
length(DrugMechDB)
x=DrugMechDB[[999]]
names(x)
x$nodes




##### drugBank: universal databases for drugs. Here, we extract drug-disease network
#rm(list=ls())
library("xml2")

path="Data/DrugBank/full_database_04Jan2023.xml"

dbdb=read_xml(path)
#kids <- xml_children(dbdb)
#memory.limit(size=1e10) #increase memory in windows
#length(kids)

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


#for (d in 1:700) res=getIndication(d)

indicationDbList=lapply(1:length(DrugBank),getIndication)
indicationDb=do.call(rbind,indicationDbList)
dim(indicationDb)
colnames(indicationDb)=c("DrugName","Description","Indication")

#save(DrugBank, indicationDb, file="DrugBank_drug_indication.RData")

########## now extract disease names from the indication
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


# ##### MESH databases
# ########### download mesh data version 2023: https://www.nlm.nih.gov/databases/download/mesh.html
# rm(list=ls())
# library("xml2")
# #load("C:/Users/truvu/OneDrive\ -\ Karolinska Institutet/Box/KI_workspace/Working/AML/Drug_repurposing/DrugBank/indicationDb_DrugBank_HPO.RData")
# 
# path="Data/desc2023.gz"
# dbdb=read_xml(path)
# 
# Mesh_xml=xml_contents(dbdb)
# length(Mesh_xml)
# # Mesh_xml=xml_contents(dbdb) #this is faster
# # length(Mesh_xml)
# # d=1
# # x=Mesh_xml[[d]]
# # x1=xml_contents(x)
# # p=6
# # x2=xml_contents(x1[p])
# # x3=xml_contents(x2[1])
# # x3
# 
# Mesh <- xml_children(dbdb) #this is slower but more structured
# length(Mesh) # 30454 cases
# d=1000
# x=as_list(Mesh[[d]])
# names(x)
# x$DescriptorUI #id
# unlist(x$DescriptorName) #name of drug/disease
# #x$ConceptList
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
# save(Mesh_des, Mesh_id, Mesh, Mesh_xml, file="Mesh.RData")



##### Human symptoms-disease network https://www.nature.com/articles/ncomms5212
### this database focuses on cancer disease, using Mesh id
library(data.table)
hsdn_net=fread(file="network_Human_symptoms–disease.txt")
hsdn_disease=fread("diseases_Human_symptoms–disease.txt")
hsdn_symptom=fread("symptoms_Human_symptoms–disease.txt")


head(hsdn_net)
head(hsdn_disease)
head(hsdn_symptom)

save(hsdn_net,hsdn_disease, hsdn_symptom, file="hsdn.RData")

####### done with all databases
###################################


##### do some matching, really big issues
rm(list=ls())

library("xml2")
library("data.table")
library("ontologyIndex")
library("ontologySimilarity")
library("ontologyPlot")


load("orphaNet.RData")
load("HPO.RData")
load("DrugBank_drug_indication_mapping.RData")
load("hsdn.RData")
load("Mesh.RData")




### need to link again to the Mesh database
path="Data/Mesh/desc2023.gz"
dbdb=read_xml(path)
Mesh_xml=xml_contents(dbdb)
Mesh <- xml_children(dbdb) #this is slower but more structured
length(Mesh) # 30454 cases

p=which(Mesh_id=="D038261")#Alexander Disease
p=which(Mesh_id=="D020739")#Brain Diseases, Metabolic, Inborn

p
x=as_list(Mesh[[p]])
unlist(x$DescriptorName) 
Mesh_des[p]

#Mesh_id[43]; Mesh_des[43] #"Abstracting and Indexing" --> meaningless
#x$TreeNumberList # each disease can belong to multiple trees --> need a map between the trees to the Mesh_id

#now match between HPO disease and MESH
hpo_Disease=toupper(unique(oMat$Disease))
Mesh_des2=toupper(Mesh_des)

hpo_mesh_map=list()
for (i in 1:length(hpo_Disease)){ #15235 drugs
  if (i %% 100 == 0) cat (" ",i)
  x1=hpo_Disease[i]
  p=grep(x1,Mesh_des2) # if x1 is included in the names: slower but can be used to check discordant names, even though only 375 matches
  #p=which(Mesh_des2 == x1) # get exact names / much faster but only 341 exist
  hpo_mesh_map[[i]]=p
}

x=lengths(hpo_mesh_map)
table(x) #only 341 diseases can be exactly matched
head(which(x==0))

hpo_Disease[1:15]


### the case that does not match : weird name, one uses Roman style, the other use Arabic style
x1=hpo_Disease[11]
x1
p=agrep(x1,Mesh_des2)
p

head(oMat[which(toupper(oMat$Disease)==x1),])

p=which(Mesh_id=="D009087")
Mesh_des2[p]



hpo_mesh_map=list()
for (i in 1:length(hpo_Disease)){ #15235 drugs
  if (i %% 100 == 0) cat (" ",i)
  x1=hpo_Disease[i]
  p=which(Mesh_des2 == x1)
  hpo_mesh_map[[i]]=p
}

### similar terms or synonyms which contain list of TermUI such as T000903687, need to check from https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/
length(x$ConceptList$Concept$TermList)



#diseases
intersect(toupper(rareX_name), toupper(hsdn_dname$`MeSH Disease Term`))
intersect(toupper(orphanet_name), toupper(hsdn_dname$`MeSH Disease Term`))





# ##################### Disease-drug association network: #https://snap.stanford.edu/biodata/datasets/10004/10004-DCh-Miner.html
# 
# miner_db=fread("DCh-Miner_miner-disease-chemical.tsv.gz")
# head(miner_db)
# 
# 



