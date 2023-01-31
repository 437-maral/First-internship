
library(ChAMP) 
library(dplyr)
library(ggplot2)
library(reshape2)
library(lumi)
library(ggpubr)
library(data.table)
library(grid)
library(gridExtra)
setwd("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/")
rawdir <- "L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw"
source('L:/!DIAB/Masoumeh Vojood/svd_new.R')
source("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/codes/dife_pca_2021_changes.R")
source("L:/!DIAB/Masoumeh Vojood/functions/PAM.R")
source('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/codes/pam_euclidean.R')



myLoad<- champ.load(directory = rawdir,
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="EPIC"
)      

save(myLoad, file= "Load.Rdata")
load('Load.Rdata')


champ.QC(Rplot=FALSE, resultsDir="./QC1/")

myNorm <- champ.norm(arraytype="EPIC",cores =1)

save(myNorm, file= "Norm.Rdata")

load("Norm.Rdata")
####

myLoad$pd$SEX='Female'
myLoad$pd$SEX=as.factor(myLoad$pd$SEX)
myLoad$pd$Sample_Name=as.factor(myLoad$pd$Sample_Name)
myLoad$pd$cluster_tueftulip=as.factor(myLoad$pd$cluster_tueftulip)
myLoad$pd$Slide=as.factor(myLoad$pd$Slide)






champ.SVD2(beta=data.frame(myNorm),pd=myLoad$pd)



myRefbase <- champ.refbase(beta = myNorm, arraytype="450K") 
champ.SVD2(beta = data.frame(myRefbase$CorrectedBeta),pd = myLoad$pd)


cellcounts_champ <- myRefbase$CellFraction
cellcounts_champ_melted <- melt(cellcounts_champ)
#levels(cellcounts_champ_melted$X2)[levels(cellcounts_champ_melted$X2)=="Gran"] <- "Neu" #label as for SeSame and Meffil 
cellcounts_champ_melted$Var1 <- as.character(cellcounts_champ_melted$Var1)
stacked_plot_cellcount_champ <- ggplot(cellcounts_champ_melted, aes(Var1, value, fill = Var2)) +
  geom_bar(stat = "identity") + 
  xlab("sample") + 
  ylab("Estimated proportion of cells")+
  ggtitle('Reference based cell type deconvolution per sample')+
  scale_x_discrete(labels= NULL)+
  scale_fill_manual(values = c("#c51b7d", "#e9a3c9", "#fde0ef", "#e6f5d0", "#a1d76a", "#4d9221"),name = "Cell type")+
  theme(text = element_text(size = 25))














new_mycombat<- champ.runCombat(beta=myRefbase$CorrectedBeta,pd=myLoad$pd,variablename=c("cluster_tueftulip"),batchname=c('Array','Slide'))
champ.SVD2(beta =data.frame(new_mycombat),pd = myLoad$pd)



save(new_mycombat, file= "new_combat.Rdata")
load("new_combat.Rdata")


###########pam clustering


#replication cohort
com=as.data.frame(new_mycombat)
pheno<-fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno.csv",skip=1,data.table=FALSE)

load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")




# filter ------------------------------------------------------------------

replication_2020<- com%>% filter(rownames(.) %in% rownames(my_combat_dmp_high_frq)) #1980
colnames(replication_2020)

replication_2020=data.frame(replication_2020)
names3 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 6)])



colnames(replication_2020)[which(colnames(replication_2020) %in% names3)]<-paste0(colnames(replication_2020)[which(colnames(replication_2020) %in% names3)],"_3")
colnames(replication_2020)[which(colnames(replication_2020) %in% names5)]<-paste0(colnames(replication_2020)[which(colnames(replication_2020) %in% names5)],"_5")
colnames(replication_2020)[which(colnames(replication_2020) %in% names6)]<-paste0(colnames(replication_2020)[which(colnames(replication_2020) %in% names6)],"_6")
colnames(replication_2020)

pam_function(input=replication_2020,class="sofya") #55.5


pam__euclidean(input=replication_2020,class="unknown3")####53 !!!!!





pam_function(input=rbind(replication_2020,pheno_short_t),class="unknown4")


pam__euclidean(input=rbind(replication_2020,pheno_short_t),class="unknown4")






















###### do pca 


cluster_rep=c(colnames(com))
cluster3_rep<- c(pheno$Sample_Name[which(pheno$cluster_tueftulip ==  '3' )])
cluster5_rep<- c(pheno$Sample_Name[which(pheno$cluster_tueftulip ==  '5')])
cluster6_rep=c(pheno$Sample_Name[which(pheno$cluster_tueftulip ==  '6')])

######
cluster_rep[which(colnames(com) %in% cluster3_rep)]<-'3'
cluster_rep[which(colnames(com) %in% cluster5_rep)]<-'5'
cluster_rep[which(colnames(com) %in% cluster6_rep)]<-'6'


PCA_rep_cluster=dife_pca(com,colorgroups = cluster_rep)



l=biplot(PCA_rep_cluster,,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Potsdam Cohort',lab = NULL)

