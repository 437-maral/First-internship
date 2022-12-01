

# libraries ----------------------------------------------------------------
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidylog)
library(janitor) #clean_names(df) or %>% clean_names()
library(ragg) #use agg_png to have quicker plots
library(report) #for session report
library(ChAMP)
require(cluster)
require(factoextra)
require(M3C)
require(tibble)
require(cvms)
require(ggthemes)
require(ggimage)
require(rsvg)

source("L:/!DIAB/Masoumeh Vojood/functions/PAM.R")
source("L:/!DIAB/Masoumeh Vojood/functions/min_max_norm.R")
source('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/codes/pam_euclidean.R')


setwd("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/")
load('Combat.Rdata')   #### combat dataset after adding phenotype

load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")
pheno<-fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno1.csv",skip=1,data.table=FALSE)

silhouette_plotter <- function(df,title="Sihouette"){
  p <- fviz_nbclust(t(df), pam, method = "silhouette")+
    theme_classic(base_size = 28)+
    ggtitle(title)+
    geom_line(aes(group = 1), color = "Purple", size = 3)
  return(p)
}




load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")



myCombat=data.frame(myCombat)
replication_2020=myCombat%>% filter(rownames(.) %in% rownames(my_combat_dmp_high_frq))##1977,48




names3 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 6)])



####colnames are just samples , so we need to rename them regarding to cluster's name
colnames(replication_2020)[which(colnames(replication_2020) %in% names3)]<-paste0(colnames(replication_2020)[which(colnames(replication_2020) %in% names3)],"_3")
colnames(replication_2020)[which(colnames(replication_2020) %in% names5)]<-paste0(colnames(replication_2020)[which(colnames(replication_2020) %in% names5)],"_5")
colnames(replication_2020)[which(colnames(replication_2020) %in% names6)]<-paste0(colnames(replication_2020)[which(colnames(replication_2020) %in% names6)],"_6")
colnames(replication_2020)




pam_function(input=replication_2020,class="sofya") #29
pam_function(input=replication_2020,class="markus")#37.72
pam_function(input=replication_2020,class="unknown1")#39
pam_function(input=replication_2020,class="unknown2")#22
pam_function(input=replication_2020,class="unknown3") #22
pam_function(input=replication_2020,class="unknown4")#47.92




pam__euclidean(input=replication_2020,class="sofya") #29
pam__euclidean(input=replication_2020,class="markus")#37.72
pam__euclidean(input=replication_2020,class="unknown1")#39
pam__euclidean(input=replication_2020,class="unknown2")#22
pam__euclidean(input=replication_2020,class="unknown3")####52 !!!!!
pam__euclidean(input=replication_2020,class="unknown4")


##adding pheno

pheno_short<-pheno %>%
  mutate(pheno_names_new=paste0("X",Sample_Name,"_",cluster_tueftulip)) %>%
  select(pheno_names_new,age,BMI,WAIST,HIP,`CRP_[mg/l]`,`HDL_[mmol/l]`,`LDL_[mmol/l]`,`Cholesterol_[mmol/l]`,`GOT(AST)_[U/l]`,`GGT_[U/l]`,`Creatinine_[µmol/l]`,MR_LFadj,`MR_SCAT/l`,`MR_VAT/l`,DI.calc,ISI.calc)


###handling missing value
pheno_short$`MR_SCAT/l`[is.na(pheno_short$`MR_SCAT/l`)] <- mean(pheno_short$`MR_SCAT/l`,na.rm = TRUE)
pheno_short$`MR_VAT/l`[is.na(pheno_short$`MR_VAT/l`)] <- mean(pheno_short$`MR_VAT/l`,na.rm = TRUE)



#####normalising
pheno_short$age<-min_max_norm(pheno_short$age)
pheno_short$BMI<-min_max_norm(pheno_short$BMI)
pheno_short$WAIST<-min_max_norm(pheno_short$WAIST)
pheno_short$HIP<-min_max_norm(pheno_short$HIP)
pheno_short$`CRP_[mg/l]`<-min_max_norm(pheno_short$`CRP_[mg/l]`)
pheno_short$`HDL_[mmol/l]`<-min_max_norm(pheno_short$`HDL_[mmol/l]`)
pheno_short$`Cholesterol_[mmol/l]`<-min_max_norm(pheno_short$`Cholesterol_[mmol/l]`)
pheno_short$`LDL_[mmol/l]`<-min_max_norm(pheno_short$`LDL_[mmol/l]`)
pheno_short$`GOT(AST)_[U/l]`<-min_max_norm(pheno_short$`GOT(AST)_[U/l]`)
pheno_short$`GGT_[U/l]`<-min_max_norm(pheno_short$`GGT_[U/l]`)
pheno_short$`Creatinine_[µmol/l]`<-min_max_norm(pheno_short$`Creatinine_[µmol/l]`)
pheno_short$MR_LFadj<-min_max_norm(pheno_short$MR_LFadj)
pheno_short$`MR_SCAT/l`<-min_max_norm(pheno_short$`MR_SCAT/l`)
pheno_short$`MR_VAT/l`<-min_max_norm(pheno_short$`MR_VAT/l`)
pheno_short$DI.calc<-min_max_norm(pheno_short$DI.calc)
pheno_short$ISI.calc<-min_max_norm(pheno_short$ISI.calc)

colnames(pheno_short)

rownames(pheno_short)<-pheno_short$pheno_names_new
pheno_short$pheno_names_new<-NULL

pheno_short_t<-t(pheno_short)


colnames(pheno_short_t)
colnames(replication_2020)
identical(colnames(pheno_short_t),colnames(replication_2020))


pam_function(input=rbind(replication_2020,pheno_short_t),class="sofya")### 54.17"



pam__euclidean(input=rbind(replication_2020,pheno_short_t),class="sofya")### 54.17"
print(silhouette_plotter(rbind(replication_2020,pheno_short_t),"2020CpGs"))




######work with intersecting region

load('Dmptest.RData')

dmp3_5=dmptest$`3_to_5`[dmptest$`3_to_5`$P.Value<0.05,]
dmp3_6=dmptest$`3_to_6`[dmptest$`3_to_6`$P.Value<0.05,]
dmp5_6=dmptest$`6_to_5`[dmptest$`6_to_5`$P.Value<0.05,]

##add Name column for cpgs
dmp3_5$Name=rownames(dmp3_5)
dmp3_6$Name=rownames(dmp3_6)
dmp5_6$Name=rownames(dmp5_6)



new=Reduce(intersect,list(dmp3_5$Name,dmp5_6$Name,dmp3_6$Name)) #1060

ident_full=myCombat[rownames(myCombat) %in% new,]

names3 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 6)])


colnames(ident_full)[which(colnames(ident_full) %in% names3)]<-paste0(colnames(ident_full)[which(colnames(ident_full) %in% names3)],"_3")
colnames(ident_full)[which(colnames(ident_full) %in% names5)]<-paste0(colnames(ident_full)[which(colnames(ident_full) %in% names5)],"_5")
colnames(ident_full)[which(colnames(ident_full) %in% names6)]<-paste0(colnames(ident_full)[which(colnames(ident_full) %in% names6)],"_6")
colnames(ident_full)


pam_function(input=ident_full,class="sofya") #93.75"
pam__euclidean(input=ident_full,class="sofya")

pam__euclidean(input=rbind(ident_full,pheno_short_t),class="sofya") #97%

##### check replication cohort with discovery cohort

inner_join(rownames_to_column(ident_full),rownames_to_column(my_combat_dmp_high_frq),by="rowname")

#> rows only in x  (1,057)
#> rows only in y  (2,017)
#> matched rows         3
#>                 =======
#> rows total           3

