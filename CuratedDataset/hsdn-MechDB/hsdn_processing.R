rm(list=ls())
library("xml2")
setwd("C:/Trang/KIProjects/hsdn/DataPreparation")

##-------------------------- Mesh database----------------------------------
path="Data/Mesh/desc2023.gz"
dbdb=read_xml(path)

Mesh_xml=xml_contents(dbdb)
length(Mesh_xml)

Mesh <- xml_children(dbdb) 
length(Mesh) # 30454 cases
d=1000
x=as_list(Mesh[[d]])
names(x)
x$DescriptorUI #id
unlist(x$DescriptorName) #name of drug/disease
#x$ConceptList

Mesh_id=sapply(Mesh_xml,function(x){
  x1=xml_contents(x)
  x2=xml_contents(x1[1])#DescriptorUI
  x3=as.character(x2)
  return(x3)
})

Mesh_des=sapply(Mesh_xml,function(x){
  x1=xml_contents(x)
  x2=xml_contents(x1[2])#DescriptorName
  x3=as.character(xml_contents(x2))
  return(x3)
})

Mesh_entryTerms=lapply(Mesh_xml,function(x){
  #  for (i in 1:length(Mesh_xml)){
  #  x=Mesh_xml[[i]]
  x1=xml_contents(x)
  k=grep("<ConceptList>", x1)
  x2=xml_contents(xml_contents(x1[k])) #extract two levels
  p=grep("<ConceptName>", x2)
  r1=as.character(xml_contents(xml_contents(x2[p]))) #level 1
  
  p=grep("<TermList>", x2)
  x3=xml_contents(xml_contents(x2[p]))
  p1=grep("<String>", x3)
  r2=as.character(xml_contents(x3[p1])) #level 2, also includes all terms in level 1
  #  }
  
  res=list(entryTerms_level1=r1, entryTerms_level2=r2)
  return(res)
})
names(Mesh_entryTerms)=Mesh_id

# save(Mesh_entryTerms, Mesh_des, Mesh_id, Mesh, Mesh_xml, file="Pre_datasets_hsdn/Mesh.RData")

load("Mesh.RData")
#combine level1 and level2 of the entry terms
Mesh_entryTerms_all=sapply(Mesh_entryTerms,function(x) unique(unlist(x)))


##-----------------------------Human symptoms-disease network (hsdn) database----------------

hsdn_net=fread("Data/network_Human_symptoms–disease.txt")
hsdn_disease=fread("Data/diseases_Human_symptoms–disease.txt")
hsdn_symptom=fread("Data/symptoms_Human_symptoms–disease.txt")


head(hsdn_net)
head(hsdn_disease)
head(hsdn_symptom)

# save(hsdn_net,hsdn_disease, hsdn_symptom, file="Pre_datasets_hsdn/hsdn.RData")

##-------------Processing hsdn and mapping with mesh id

load("Pre_datasets_hsdn/hsdn.RData")
load("Pre_datasets_hsdn/Mesh.RData")

Mesh_id2=Mesh_id
names(Mesh_id2)=Mesh_des

Mesh_des2=Mesh_des
names(Mesh_des2)=Mesh_id

hsdn_symptom$Mesh_id=Mesh_id2[hsdn_symptom$`MeSH Symptom Term`]

m=hsdn_symptom$`MeSH Symptom Term` %in% Mesh_des
table(m)
m=which(!m)
notFoundsymptom=hsdn_symptom[m,]
dim(notFoundsymptom)

for (i in 1:nrow(notFoundsymptom)){
  s=notFoundsymptom$`MeSH Symptom Term`[i]
  y=sapply(Mesh_entryTerms_all, function(x) which(s %in% x))
  j=which(lengths(y)>0)
  if (length(j)>0) notFoundsymptom$Mesh_id[i]=names(Mesh_entryTerms_all[j])
}

hsdn_symptom$Mesh_id[m]=notFoundsymptom$Mesh_id
hsdn_symptom$Mesh_des=Mesh_des2[hsdn_symptom$Mesh_id]

#found in the entry terms, replace by the representative ones
#notFoundsymptom$Replacement=c("Hypoxia","Intellectual Disability","Sleep Wake Disorders", "Learning Disabilities","IgA Vasculitis","Alien Limb Phenomenon","Brief, Resolved, Unexplained Event") 


hsdn_disease$Mesh_id=Mesh_id2[hsdn_disease$`MeSH Disease Term`]

m=hsdn_disease$`MeSH Disease Term` %in% Mesh_des
table(m)
m=which(!m)
notFounddisease=hsdn_disease[m,]
dim(notFounddisease)
#notFounddisease #some symptoms are the names of diseases as well
for (i in 1:nrow(notFounddisease)){
  s=notFounddisease$`MeSH Disease Term`[i]
  y=sapply(Mesh_entryTerms_all, function(x) which(s %in% x))
  j=which(lengths(y)>0)
  cat(" ",length(j))
  if (length(j)>0)
    notFounddisease$Mesh_id[i]=names(Mesh_entryTerms_all[j])
}
hsdn_disease$Mesh_id[m]=notFounddisease$Mesh_id
hsdn_disease$Mesh_des=Mesh_des2[hsdn_disease$Mesh_id]

p=which(is.na(hsdn_disease$Mesh_id)) #some dieases do not match with the current database
p
hsdn_disease[p,]
hsdn_disease$Mesh_id[3200]= "D007863"
hsdn_disease$Mesh_id[4299]= "D000080443"
hsdn_disease$Mesh_id[2249]= "D000083682" #can be drug overdose or Opiate Overdose, I selected Opiate Overdose because it is more specific

hsdn_disease$Mesh_des=Mesh_des2[hsdn_disease$Mesh_id]


dim(hsdn_symptom)
dim(hsdn_disease)
dim(hsdn_net)

# save(hsdn_net,hsdn_disease, hsdn_symptom, file="Pre_datasets_hsdn/hsdn_withMeshId.RData")
