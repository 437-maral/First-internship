
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


myLoad <- champ.load(directory = rawdir,
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

myLoad$pd$SEX='Female'
myLoad$pd$SEX=as.factor(myLoad$pd$SEX)
myLoad$pd$Sample_Name=as.factor(myLoad$pd$Sample_Name)
myLoad$pd$cluster_tueftulip=as.factor(myLoad$pd$cluster_tueftulip)
myLoad$pd$Slide=as.factor(myLoad$pd$Slide)

#norm
myNorm <- champ.norm(arraytype="EPIC",cores =1)

save(myNorm, file= "Norm.Rdata")


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



myCombat <- champ.runCombat(beta=myRefbase$CorrectedBeta,pd=myLoad$pd,variablename=c("cluster_tueftulip"),batchname=c('Slide','Array'))
champ.SVD2(beta = data.frame(myCombat),pd = myLoad$pd)

save(myCombat,file='Combat.Rdata')

load('Load.Rdata')


dmptest <- champ.DMP(beta = myCombat,
                      pheno = as.character(myLoad$pd$cluster_tueftulip),
                      compare.group = NULL,
                      adjPVal = 1,
                      adjust.method = "BH",
                      arraytype = "EPIC")

save(dmptest,file = "Dmptest.Rdata")


load('Dmptest.Rdata')

