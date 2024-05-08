#### doing pca , to see outlires
#-------library
library(tidyverse)
library(broom)
library(data.table)
library(PCAtools)
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidylog)

#---------function
source("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/codes/dife_pca_2021_changes.R")

#-------------file
setwd('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort')
###phnotypes of discovery Cohort
pheno1=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort//pheno.csv',skip = 1,data.table = FALSE)
pheno1$Sample_Name <- as.character(pheno1$Sample_Name)#like last time 
pheno1$cluster_tueftulip <- as.factor(pheno1$cluster_tueftulip)

#phenotypes für replication cohort
pheno2<-fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno1.csv",skip = 1,data.table = FALSE)
pheno2$Sample_Name <- as.character(pheno2$Sample_Name)#like last time 
pheno2$cluster_tueftulip <- as.factor(pheno2$cluster_tueftulip)

#######
load('Combat.Rdata')
mycombat_rep=myCombat    


load('femonly_3rm_ref_plate.Rdata')
mycombat_disc=myCombat
#####




###############repliecation Cohort
########## colors baed on batch effects



ID_rep=c(colnames(mycombat_rep))
####slide
#replication Cohort
groupA <- c(pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026050025)])
groupB <- c(pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026050032)])
groupC<- c(pheno2$Sample_Name[which(pheno2$Sentrix_ID == 206026060040)])
groupD<- c(pheno2$Sample_Name[which(pheno2$Sentrix_ID == 206026060052)])
groupE<- c(pheno2$Sample_Name[which(pheno2$Sentrix_ID == 206026060053)])
groupF<- c(pheno2$Sample_Name[which(pheno2$Sentrix_ID == 206026060054)])
###
ID_rep[which(colnames(mycombat_rep) %in% groupA)]<-'GroupA'
ID_rep[which(colnames(mycombat_rep) %in% groupB)]<-'GroupB'
ID_rep[which(colnames(mycombat_rep) %in% groupC)]<-'GroupC'
ID_rep[which(colnames(mycombat_rep) %in% groupD)]<-'GroupD'
ID_rep[which(colnames(mycombat_rep) %in% groupE)]<-'GroupE'
ID_rep[which(colnames(mycombat_rep) %in% groupF)]<-'GroupF'


#discovery Cohort
disc_ID=c(colnames(mycombat_disc))
slideA <- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624640042)])
slideB <- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624640068)])
slideC<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624640048)])
slideD<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624640064)])
slideE<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624640067)])
slideF<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624880101)])
slideG<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624640031)])
slideH<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624870177)])
slideI<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624860150)])
slideM<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624640041)])
slideN<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624640070)])
slideV<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624870084)])
slideZ<<- c(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624880134)])


disc_ID[which(colnames(mycombat_disc) %in% slideA)]<-'GroupA'
disc_ID[which(colnames(mycombat_disc) %in% slideB)]<-'GroupB'
disc_ID[which(colnames(mycombat_disc) %in% slideC)]<-'GroupC'
disc_ID[which(colnames(mycombat_disc) %in% slideD)]<-'GroupD'
disc_ID[which(colnames(mycombat_disc) %in% slideE)]<-'GroupE'
disc_ID[which(colnames(mycombat_disc) %in% slideF)]<-'GroupF'

disc_ID[which(colnames(mycombat_disc) %in% slideG)]<-'GroupG'
disc_ID[which(colnames(mycombat_disc) %in% slideH)]<-'GroupH'
disc_ID[which(colnames(mycombat_disc) %in% slideI)]<-'GroupI'
disc_ID[which(colnames(mycombat_disc) %in% slideM)]<-'GroupM'
disc_ID[which(colnames(mycombat_disc) %in% slideN)]<-'GroupN'
disc_ID[which(colnames(mycombat_disc) %in% slideV)]<-'GroupV'
disc_ID[which(colnames(mycombat_disc) %in% slideZ)]<-'GroupZ'


#####I changed the pca function , to  get  only pca result 
PCA_rep_id=dife_pca(mycombat_rep,colorgroups = ID_rep)

PCA_disc_id=dife_pca(mycombat_disc,colorgroups = disc_ID)

pdf('PCA_slide.pdf')
print(biplot(PCA_rep_id,lab=colnames(mycombat_rep),colby ='group',hline = 0, vline = 0,legendPosition = 'right'))
print(biplot(PCA_disc_id,lab=colnames(mycombat_disc),colby ='group',hline = 0, vline = 0,legendPosition = 'right'))
dev.off()



######array
array_rep=c(colnames(mycombat_rep))
arrayA <- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R01C01'  )])
arrayB<- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R02C01')])
ArrayC<- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position == 'R03C01')])
ArrayD<- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R04C01')])
arrayE <- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R05C01'  )])
arrayF<- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R06C01')])
ArrayG<- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position == 'R07C01')])
ArrayH<- c(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R08C01')])

array_rep[which(colnames(mycombat_rep) %in% arrayA)]<-'GroupA'
array_rep[which(colnames(mycombat_rep) %in% arrayB)]<-'GroupB'
array_rep[which(colnames(mycombat_rep) %in% ArrayC)]<-'GroupC'
array_rep[which(colnames(mycombat_rep) %in% ArrayD)]<-'GroupD'
array_rep[which(colnames(mycombat_rep) %in% arrayE)]<-'GroupE'
array_rep[which(colnames(mycombat_rep) %in% arrayF)]<-'GroupF'
array_rep[which(colnames(mycombat_rep) %in% ArrayG)]<-'GroupG'
array_rep[which(colnames(mycombat_rep) %in% ArrayH)]<-'GroupH'

###discovery Cohort


array_disc=c(colnames(mycombat_disc))
Array1 <- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R01C01'  )])
Array2<- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R02C01')])
Array3<- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position == 'R03C01')])
Array4<- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R04C01')])
Array5 <- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R05C01'  )])
Array6<- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R06C01')])
Array7<- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position == 'R07C01')])
Array8<- c(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R08C01')])

array_disc[which(colnames(mycombat_disc) %in% Array1 )]<-'GroupA'
array_disc[which(colnames(mycombat_disc) %in% Array2)]<-'GroupB'
array_disc[which(colnames(mycombat_disc) %in% Array3 )]<-'GroupC'
array_disc[which(colnames(mycombat_disc) %in% Array4)]<-'GroupD'
array_disc[which(colnames(mycombat_disc) %in% Array5 )]<-'GroupE'
array_disc[which(colnames(mycombat_disc) %in% Array6 )]<-'GroupF'
array_disc[which(colnames(mycombat_disc) %in% Array7)]<-'GroupG'
array_disc[which(colnames(mycombat_disc) %in% Array8 )]<-'GroupH'


PCA_rep_array=dife_pca(mycombat_rep,colorgroups = array_rep)

PCA_disc_array=dife_pca(mycombat_disc,colorgroups = array_disc)

pdf('PCA_array.pdf')
print(biplot(PCA_rep_array,lab=colnames(mycombat_rep),colby ='group',hline = 0, vline = 0,legendPosition = 'right'))
print(biplot(PCA_disc_array,lab=colnames(mycombat_disc),colby ='group',hline = 0, vline = 0,legendPosition = 'right'))
dev.off()




####### cluster 


##replication Cohort
cluster_rep=c(colnames(mycombat_rep))
cluster3_rep<- c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '3' )])
cluster5_rep<- c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '5')])
cluster6_rep=c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '6')])

######
cluster_rep[which(colnames(mycombat_rep) %in% cluster3_rep)]<-'3'
cluster_rep[which(colnames(mycombat_rep) %in% cluster5_rep)]<-'5'
cluster_rep[which(colnames(mycombat_rep) %in% cluster6_rep)]<-'6'



####discovery Cohort

cluster_disc=c(colnames(mycombat_disc))
cluster3_disc<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '3' )])
cluster5_disc<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '5')])
cluster6_disc<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '6')])

######
cluster_disc[which(colnames(mycombat_disc) %in% cluster3_disc)]<-'3'
cluster_disc[which(colnames(mycombat_disc) %in% cluster5_disc)]<-'5'
cluster_disc[which(colnames(mycombat_disc) %in% cluster6_disc)]<-'6'


PCA_rep_cluster=dife_pca(mycombat_rep,colorgroups = cluster_rep)

PCA_disc_cluster=dife_pca(mycombat_disc,colorgroups = cluster_disc)


l=biplot(PCA_rep_cluster,,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Potsdam Cohort',lab = NULL)
f=biplot(PCA_disc_cluster,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Tübingen Cohort',lab = NULL,xlim = c(-10,10),ylim=c(-10,10))

library(ggpubr)
ggarrange(l,f)


#-----------------------------------------

###### look at merging of two cohorts

mycombat_disc=data.frame(mycombat_disc)
mycombat_rep=data.frame(mycombat_rep)
####
join_twochorts=merge(rownames_to_column(mycombat_disc),rownames_to_column(mycombat_rep))
###
rownames(join_twochorts)=join_twochorts$rowname
join_twochorts=select(join_twochorts,-rowname)

####color based on cluster 


two_cohorts=c(colnames(join_twochorts))

'discovery Cohort'
names3_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 3)])
names5_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 5)])
names6_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 6)])

#replication cohort
names3_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 3)])
names5_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 5)])
names6_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 6)])

######

two_cohorts[which(colnames(join_twochorts) %in% names3_disc)]<-'3'
two_cohorts[which(colnames(join_twochorts) %in% names5_disc)]<-'5'
two_cohorts[which(colnames(join_twochorts) %in% names6_disc)]<-'6'



two_cohorts[which(colnames(join_twochorts) %in% names3_rep)]<-'3'
two_cohorts[which(colnames(join_twochorts) %in% names5_rep)]<-'5'
two_cohorts[which(colnames(join_twochorts) %in% names6_rep)]<-'6'




metadata_cohorts<-data.frame(row.names=colnames(join_twochorts))
metadata_cohorts$group<-two_cohorts



pca_cohort<- pca(join_twochorts,metadata = metadata_cohorts)


biplot(pca_cohort,lab=NULL,colby ='group',hline = 0, vline = 0,legendPosition = 'right')


ggsave('cohort.png')



####regarding slide 
slide_twocohort=c(colnames(join_twochorts))


slide1_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026050025)]))
slide2_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026050032)])
slide3_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_ID == 206026060040)]))
slide4_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026060052)]))
slide5_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026060053)])
slide6_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_ID ==  206026060054)]))


slide1_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==   205624640042)]))
slide2_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_ID ==  205624640068)])
slide3_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==  205624640048)]))
slide4_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==   205624640064 )]))
slide5_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_ID ==  205624640067)])
slide6_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==  205624880101)]))
slide7_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==   205624640031)]))
slide8_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_ID ==  205624870177)])
slide9_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624860150)]))
slide10_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==  205624640041 )]))
slide11_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624640070)])
slide12_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID ==205624870084)]))
slide13_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_ID == 205624880134)]))



slide_twocohort[which(colnames(join_twochorts) %in% slide1_rep | colnames(join_twochorts) %in% slide1_disc )]<-'GroupA'
slide_twocohort[which(colnames(join_twochorts) %in% slide2_rep | colnames(join_twochorts) %in% slide2_disc )]<-'GroupB'
slide_twocohort[which(colnames(join_twochorts) %in% slide3_rep | colnames(join_twochorts) %in% slide3_disc )]<-'GroupC'
slide_twocohort[which(colnames(join_twochorts) %in% slide4_rep | colnames(join_twochorts) %in% slide4_disc )]<-'GroupD'
slide_twocohort[which(colnames(join_twochorts) %in% slide5_rep | colnames(join_twochorts) %in% slide5_disc )]<-'GroupE'
slide_twocohort[which(colnames(join_twochorts) %in% slide6_rep | colnames(join_twochorts) %in% slide6_disc )]<-'GroupF'

slide_twocohort[which(colnames(join_twochorts) %in% slide7_disc )]<-'GroupG'
slide_twocohort[which(colnames(join_twochorts) %in% slide8_disc )]<-'GroupH'
slide_twocohort[which(colnames(join_twochorts) %in% slide9_disc )]<-'GroupI'
slide_twocohort[which(colnames(join_twochorts) %in% slide10_disc )]<-'GroupJ'
slide_twocohort[which(colnames(join_twochorts) %in% slide11_disc )]<-'GroupK'
slide_twocohort[which(colnames(join_twochorts) %in% slide12_disc )]<-'GroupL'
slide_twocohort[which(colnames(join_twochorts) %in% slide13_disc )]<-'GroupM'



metadata_slide<-data.frame(row.names=colnames(join_twochorts))
metadata_slide$group<-slide_twocohort


pca_slide<- pca(join_twochorts,metadata = metadata_slide)


####array 

array_twocohort=c(colnames(join_twochorts))

array1_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R01C01')]))
array2_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R02C01')])
array3_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_Position == 'R03C01')]))
array4_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R04C01')]))
array5_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R05C01')])
array6_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_Position == 'R06C01')]))
array7_rep<- paste0("X",pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R07C01')])
array8_rep<- paste0("X",(pheno2$Sample_Name[which(pheno2$Sentrix_Position ==  'R08C01')]))



array1_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R01C01')]))
array2_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R02C01')])
array3_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_Position == 'R03C01')]))
array4_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R04C01')]))
array5_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R05C01')])
array6_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_Position == 'R06C01')]))
array7_disc<- paste0("X",pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R07C01')])
array8_disc<- paste0("X",(pheno1$Sample_Name[which(pheno1$Sentrix_Position ==  'R08C01')]))



array_twocohort[which(colnames(join_twochorts) %in% array1_rep )]<-'GroupA'
array_twocohort[which(colnames(join_twochorts) %in% array2_rep )]<-'GroupB'
array_twocohort[which(colnames(join_twochorts) %in% array3_rep )]<-'Groupc'
array_twocohort[which(colnames(join_twochorts) %in% array4_rep )]<-'GroupD'
array_twocohort[which(colnames(join_twochorts) %in% array5_rep )]<-'GroupE'
array_twocohort[which(colnames(join_twochorts) %in% array6_rep )]<-'GroupF'
array_twocohort[which(colnames(join_twochorts) %in% array7_rep )]<-'GroupG'
array_twocohort[which(colnames(join_twochorts) %in% array8_rep )]<-'GroupH'


array_twocohort[which(colnames(join_twochorts) %in% array1_disc )]<-'GroupA'
array_twocohort[which(colnames(join_twochorts) %in% array2_disc )]<-'GroupB'
array_twocohort[which(colnames(join_twochorts) %in% array3_disc)]<-'Groupc'
array_twocohort[which(colnames(join_twochorts) %in% array4_disc)]<-'GroupD'
array_twocohort[which(colnames(join_twochorts) %in% array5_disc )]<-'GroupE'
array_twocohort[which(colnames(join_twochorts) %in% array6_disc)]<-'GroupF'
array_twocohort[which(colnames(join_twochorts) %in% array7_disc)]<-'GroupG'
array_twocohort[which(colnames(join_twochorts) %in% array8_disc )]<-'GroupH'



metadata_array<-data.frame(row.names=colnames(join_twochorts))
metadata_array$group<-array_twocohort

pca_array<- pca(join_twochorts,metadata = metadata_array)



pdf('PCA_cohort.pdf')
print(biplot(pca_,lab=NULL,colby ='group',hline = 0, vline = 0,legendPosition = 'right')
dev.off()



