
###in This script , I used pam_clustering on dmps from charloot_ling method 
##I got high accuracy , so I tried with  Tübingen as well


###libraries
library(cluster)
library(data.table)
library(ggVennDiagram)
##directory

setwd("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/Charlotte_Ling_method")


##function&data
source("L:/!DIAB/Masoumeh Vojood/functions/PAM.R")
source("L:/!DIAB/Masoumeh Vojood/functions/min_max_norm.R")

load('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/two_cohorts/Charlotte_Ling_method/potsdam_dmp.Rdata')
load('new_combat.Rdata')
pheno_tub=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno.csv',skip = 1,data.table = FALSE)


dmp3=dmp_uni_3$`6_to_3`[dmp_uni_3$`6_to_3`$P.Value< 0.05, ]
dmp5=dmp_uni_5$`6_to_5`[dmp_uni_5$`6_to_5`$P.Value< 0.05, ]
dmp6=dmp_uni_6$`6_to_3`[dmp_uni_6$`6_to_3`$P.Value< 0.05, ]

dmp3_adj1=dmp3[dmp3$adj.P.Val<0.05,]#5,503

dmp5_adj1=dmp5[dmp5$adj.P.Val<0.05,]#791

dmp6_adj1=dmp6[dmp6$adj.P.Val<0.05,]#2,198

dmp3_adj1$Name=rownames(dmp3_adj1)
dmp5_adj1$Name=rownames(dmp5_adj1)
dmp6_adj1$Name=rownames(dmp6_adj1)



adj1=list('3'=dmp3_adj1$Name,'5'=dmp5_adj1$Name,'6'=dmp6_adj1$Name)


ggVennDiagram(adj1,label_alpha = 0,set_size=8,label_size=10,label="count") +  ggplot2::scale_fill_gradient(low="white",high = "white")+  theme(legend.position="none")


##all_dmps
all_cpgs=unique(c(dmp3_adj1$Name,dmp5_adj1$Name,dmp6_adj1$Name))

mycombat_alldmps=new_mycombat[rownames(new_mycombat) %in% all_cpgs,]# 8,248 



names3 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 3)]
names5 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 5)]
names6 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 6)]


colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names3)]<-paste0(colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names3)],"_3")
colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names5)]<-paste0(colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names5)],"_5")
colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names6)]<-paste0(colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names6)],"_6")
colnames(mycombat_alldmps)



pam_function(input=mycombat_alldmps,class="unknown4")#97%"



dmps=c(dmp3_adj1$Name,dmp5_adj1$Name,dmp6_adj1$Name)
duplicated  <- dmps[duplicated(dmps)]



mycombat_intersecting=new_mycombat[rownames(new_mycombat) %in% duplicated ,]




names3 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 3)]
names5 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 5)]
names6 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 6)]

colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names3)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names3)],"_3")
colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names5)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names5)],"_5")
colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names6)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names6)],"_6")
colnames(mycombat_intersecting)



pam_function(input=mycombat_intersecting,class='unknown4')#97 %



duplicated  <- dmps[duplicated(dmps)] #22
specific_dmps <- subset(dmps, !dmps %in% duplicated) #65k


mycombat_spec=new_mycombat[rownames(new_mycombat) %in% specific_dmps,]#  65150 





names3 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 3)]
names5 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 5)]
names6 <- pheno_tub$Sample_Name[which(pheno_tub$cluster_tueftulip == 6)]


colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names3)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names3)],"_3")
colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names5)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names5)],"_5")
colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names6)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names6)],"_6")
colnames(mycombat_spec)



pam_function(input=mycombat_spec,class="unknown2")###7
####phenotypes




####I would like to take these cpgs and try on tübingen_ cohort

rm(dmp3,dmp5,dmp6)

rm(dmp_uni_3,dmp_uni_5,dmp_uni_6)


###work with  dmp_adjusted


rm(new_mycombat)
#####phenotypes
library(dplyr)

load()
pheno_tubingen=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/pheno.csv',skip = 1,data.table = FALSE)


load('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/Charlotte_Ling_method/Full_data_for_ling_method.Rdata')

pheno_tub=fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/pheno.csv",skip=1,data.table=FALSE)



phenotypes_cleaned <- pheno_tubingen %>%
  janitor::clean_names(
    case = "snake",
    replace = c("µ" = "micro", `'` = "", `"` = "", `%` = "_percent_", `#` = "_number_")
  ) %>% na_if(.,"n.a.") %>% #replace the n.a. with real NAs 
  mutate(cluster_whitehall = as.numeric(cluster_whitehall)) %>%
  mutate(cluster_tueftulip = as.numeric(cluster_tueftulip)) %>%
  mutate(waist = as.numeric(waist)) %>%
  mutate(hip = as.numeric(hip)) %>%
  mutate(mr_tat_l = as.numeric(mr_tat_l)) %>%
  mutate(glucose_000 = as.numeric(glucose_000)) %>%
  mutate(glucose_120 = as.numeric(glucose_120)) %>%
  mutate(glucose_auc = as.numeric(glucose_auc)) %>%
  mutate(cholesterol_mmol_l = as.numeric(cholesterol_mmol_l)) %>%
  mutate(ldl_mmol_l = as.numeric(ldl_mmol_l)) %>%
  mutate(triglycerides_mmol_l = as.numeric(triglycerides_mmol_l)) %>%
  mutate(got_ast_u_l = as.numeric(got_ast_u_l)) %>%
  mutate(gpt_alt_u_l = as.numeric(gpt_alt_u_l)) %>%
  mutate(ggt_u_l = as.numeric(ggt_u_l)) %>%
  mutate(creatinine_micromol_l = as.numeric(creatinine_micromol_l)) %>%
  mutate(igi_calc = as.numeric(igi_calc)) %>%
  dplyr::select(1:28,prs,-sex,-study_date)

easy_measure_phenos <- phenotypes_cleaned %>%
  dplyr::select("sample_name",
                "bmi","waist","hip","crp_mg_l","glucose_000","cholesterol_mmol_l",
                "hdl_mmol_l","ldl_mmol_l","triglycerides_mmol_l","got_ast_u_l","gpt_alt_u_l","ggt_u_l","creatinine_micromol_l"
  )

###mycombat = is tübingen cohort


pheno_tubingen=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/pheno.csv',skip = 1,data.table = FALSE)
tub_alldmps=myCombat[rownames(myCombat) %in% all_cpgs,] #8185 72


names3_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 3)]
names5_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 5)]
names6_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 6)]

colnames(tub_alldmps)[which(colnames(tub_alldmps) %in% names3_tub)]<-paste0(colnames(tub_alldmps)[which(colnames(tub_alldmps) %in% names3_tub)],"_3")
colnames(tub_alldmps)[which(colnames(tub_alldmps) %in% names5_tub)]<-paste0(colnames(tub_alldmps)[which(colnames(tub_alldmps) %in% names5_tub)],"_5")
colnames(tub_alldmps)[which(colnames(tub_alldmps) %in% names6_tub)]<-paste0(colnames(tub_alldmps)[which(colnames(tub_alldmps) %in% names6_tub)],"_6")
colnames(tub_alldmps)



pam_function(input=tub_alldmps,class='markus')#34%
pam_function(input=tub_alldmps,class='sofya')#26%31%
pam_function(input=tub_alldmps,class='unknown1')#31 %
pam_function(input=tub_alldmps,class='unknown2')#33 %
pam_function(input=tub_alldmps,class='unknown3')#26 %
pam_function(input=tub_alldmps,class='unknown4')#47 %

####adding phenotypes 



pheno_preper <- function (choosen_pheno,names3,names5,names6){
  #we only can take phenotypes where we have data for ALL samples, kick the rest
  choosen_pheno <- choosen_pheno[rowSums(is.na(choosen_pheno)) == 0,]
  rownames(choosen_pheno) <- choosen_pheno$sample_name
  #normalise with minmax normalisation
  
  cluster_choosen_pheno_norm <-  as.data.frame(lapply(choosen_pheno[2:ncol(choosen_pheno)], min_max_norm))
  cluster_choosen_pheno_norm_t <- t(cluster_choosen_pheno_norm)
  colnames(cluster_choosen_pheno_norm_t) <- choosen_pheno$sample_name
  colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names3] <- paste0(colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names3] ,"_3")
  colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names5] <- paste0(colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names5] ,"_5")
  colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names6] <- paste0(colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names6] ,"_6")
  
  return(cluster_choosen_pheno_norm_t)
}

pheno_easy<- easy_measure_phenos[which(easy_measure_phenos$sample_name %in% gsub("_.*","",colnames(tub_alldmps))),]



cluster_easy_pheno<- pheno_preper(pheno_easy,names3_tub,names5_tub,names6_tub)


pam_function(input=rbind(cluster_easy_pheno,tub_alldmps),class='markus')
pam_function(input=rbind(cluster_easy_pheno,tub_alldmps),class='sofya')
pam_function(input=rbind(cluster_easy_pheno,tub_alldmps),class='unknown1')
pam_function(input=rbind(cluster_easy_pheno,tub_alldmps),class='unknown2')
pam_function(input=rbind(cluster_easy_pheno,tub_alldmps),class='unknown3')
pam_function(input=rbind(cluster_easy_pheno,tub_alldmps),class='unknown4')


###intersecting region


tub_intersecting=myCombat[rownames(myCombat) %in% duplicated,]


names3_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 3)]
names5_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 5)]
names6_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 6)]

colnames(tub_intersecting)[which(colnames(tub_intersecting) %in% names3_tub)]<-paste0(colnames(tub_intersecting)[which(colnames(tub_intersecting) %in% names3_tub)],"_3")
colnames(tub_intersecting)[which(colnames(tub_intersecting) %in% names5_tub)]<-paste0(colnames(tub_intersecting)[which(colnames(tub_intersecting) %in% names5_tub)],"_5")
colnames(tub_intersecting)[which(colnames(tub_intersecting) %in% names6_tub)]<-paste0(colnames(tub_intersecting)[which(colnames(tub_intersecting) %in% names6_tub)],"_6")
colnames(tub_intersecting)



pam_function(input=tub_intersecting,class='markus')#22%
pam_function(input=tub_intersecting,class='sofya')#33%
pam_function(input=tub_intersecting,class='unknown1')#38 %
pam_function(input=tub_intersecting,class='unknown2')#38 %
pam_function(input=tub_intersecting,class='unknown3')#22 %
pam_function(input=tub_intersecting,class='unknown4')#44 %244



####

pheno_easy_int<- easy_measure_phenos[which(easy_measure_phenos$sample_name %in% gsub("_.*","",colnames(tub_intersecting))),]



cluster_easy_pheno_int<- pheno_preper(pheno_easy_int,names3_tub,names5_tub,names6_tub)


pam_function(input=rbind(cluster_easy_pheno_int,tub_intersecting),class='markus')
pam_function(input=rbind(cluster_easy_pheno_int,tub_intersecting),class='sofya')
pam_function(input=rbind(cluster_easy_pheno_int,tub_intersecting),class='unknown1')
pam_function(input=rbind(cluster_easy_pheno_int,tub_intersecting),class='unknown2')
pam_function(input=rbind(cluster_easy_pheno_int,tub_intersecting),class='unknown3')
pam_function(input=rbind(cluster_easy_pheno_int,tub_intersecting),class='unknown4')



####specefic

duplicated  <- dmps[duplicated(dmps)] #22
specific_dmps <- subset(dmps, !dmps %in% duplicated) #65k


tub_spec=myCombat[rownames(myCombat) %in% specific_dmps,]#  7941




names3_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 3)]
names5_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 5)]
names6_tub<- pheno_tubingen$Sample_Name[which(pheno_tubingen$cluster_tueftulip == 6)]


colnames(tub_spec)[which(colnames(tub_spec) %in% names3_tub)]<-paste0(colnames(tub_spec)[which(colnames(tub_spec) %in% names3_tub)],"_3")
colnames(tub_spec)[which(colnames(tub_spec) %in% names5_tub)]<-paste0(colnames(tub_spec)[which(colnames(tub_spec) %in% names5_tub)],"_5")
colnames(tub_spec)[which(colnames(tub_spec) %in% names6_tub)]<-paste0(colnames(tub_spec)[which(colnames(tub_spec) %in% names6_tub)],"_6")
colnames(tub_spec)



pam_function(input=tub_spec,class='markus')#34%
pam_function(input=tub_spec,class='sofya')#26
pam_function(input=tub_spec,class='unknown1')#31 %
pam_function(input=tub_spec,class='unknown2')#33 %
pam_function(input=tub_spec,class='unknown3')#26 %
pam_function(input=tub_spec,class='unknown4')#47 %



pheno_easy_sp<- easy_measure_phenos[which(easy_measure_phenos$sample_name %in% gsub("_.*","",colnames(tub_spec))),]



cluster_easy_pheno_sp<- pheno_preper(pheno_easy_sp,names3_tub,names5_tub,names6_tub)


pam_function(input=rbind(cluster_easy_pheno_sp,tub_spec),class='markus')
pam_function(input=rbind(cluster_easy_pheno_sp,tub_spec),class='sofya')
pam_function(input=rbind(cluster_easy_pheno_sp,tub_spec),class='unknown1')
pam_function(input=rbind(cluster_easy_pheno_sp,tub_spec),class='unknown2')
pam_function(input=rbind(cluster_easy_pheno_sp,tub_spec),class='unknown3')
pam_function(input=rbind(cluster_easy_pheno_sp,tub_spec),class='unknown4')





