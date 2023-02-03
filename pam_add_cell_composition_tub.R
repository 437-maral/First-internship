##adding cell_composition as phenotypes

library(data.table)
library(cluster)
library(tibble)
library(tidyr)
library(gridExtra)
library(factoextra)
library(ggrepel)
library(cluster)
library(dplyr)
library(stringr)
library(openxlsx)
library(cvms)
library(ggVennDiagram)


####
main_folder<-"L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/Clusterings/"
####Loading data
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/DMPs_combat.RData")
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction/DMPs_without_celltype_correction.RData")

#cell_composition---------------------------------



cell_com_potsdam=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/cellcomposition.tsv',data.table=FALSE)






source("L:/!DIAB/Masoumeh Vojood/functions/min_max_norm.R")
source("L:/!DIAB/Masoumeh Vojood/functions/PAM_2version.R")

#FUNCTION-----------------------------

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

###potsdam
dmp_3_6_Potsdam_0.05<- DMP_potsdam[["6_to_3"]] %>% filter(P.Value<0.05)
dmp_3_5__Potsdam_0.05<- DMP_potsdam[["5_to_3"]] %>% filter(P.Value<0.05)
dmp_5_6_Potsdam_0.05<- DMP_potsdam[["6_to_5"]] %>% filter(P.Value<0.05)


dmp_3_6_Potsdam_0.05$Name=rownames(dmp_3_6_Potsdam_0.05)
dmp_3_5__Potsdam_0.05$Name=rownames(dmp_3_5__Potsdam_0.05)
dmp_5_6_Potsdam_0.05$Name=rownames(dmp_5_6_Potsdam_0.05)




#####tubingen
dmp_3_6_tubingen<- DMP[["3_to_6"]] %>% filter(P.Value<0.05)
dmp_3_5_tubingen<- DMP[["3_to_5"]] %>% filter(P.Value<0.05)
dmp_5_6_tubingen<- DMP[["6_to_5"]] %>% filter(P.Value<0.05)


dmp_3_6_tubingen$Name=rownames(dmp_3_6_tubingen)
dmp_3_5_tubingen$Name=rownames(dmp_3_5_tubingen)
dmp_5_6_tubingen$Name=rownames(dmp_5_6_tubingen)



#### do overlapping
##0.05

ncpgs_3_5_0.05<- dmp_3_5__Potsdam_0.05 %>%
  filter(Name %in% dmp_3_5_tubingen$Name)

ncpgs_3_6_0.05<-dmp_3_6_Potsdam_0.05 %>%
  filter(Name %in% dmp_3_6_tubingen$Name)

ncpgs_5_6_0.05<-dmp_5_6_Potsdam_0.05 %>%
  filter(Name %in% dmp_5_6_tubingen$Name)

plt=ggVennDiagram(list("3_5"=ncpgs_3_5_0.05$Name,"3_6"=ncpgs_3_6_0.05$Name,"5_6"=ncpgs_5_6_0.05$Name),label_alpha = 0,set_size=3,label_size=3) +
  ggplot2::scale_fill_gradient(low="white",high = "white")+ggtitle('Pvalue_0.05')

only_35_36<-plt[['plot_env']][["data"]]@region[["item"]][[4]]
only_35_56<-plt[["plot_env"]][["data"]]@region[["item"]][[5]]
only_36_56<-plt[["plot_env"]][["data"]]@region[["item"]][[6]]


intersecting=c(only_35_36,only_35_56,only_36_56)#701


#tub---------------------------------------
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction/myNorm.RData")


cpgs<-myNorm

# pheno --------------------------------------------------------------------
pheno_tub=fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/pheno.csv",skip=1,data.table=FALSE)




#####with phenotypes 
###cleaning 

phenotypes_cleaned <- pheno_tub %>%
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






names3 <- phenotypes_cleaned$sample_name[which(phenotypes_cleaned$cluster_tueftulip == 3)]
names5 <- phenotypes_cleaned$sample_name[which(phenotypes_cleaned$cluster_tueftulip == 5)]
names6 <- phenotypes_cleaned$sample_name[which(phenotypes_cleaned$cluster_tueftulip == 6)]


colnames(cpgs)[which(colnames(cpgs) %in% names3)]<-paste0(colnames(cpgs)[which(colnames(cpgs) %in% names3)],"_3")
colnames(cpgs)[which(colnames(cpgs) %in% names5)]<-paste0(colnames(cpgs)[which(colnames(cpgs) %in% names5)],"_5")
colnames(cpgs)[which(colnames(cpgs) %in% names6)]<-paste0(colnames(cpgs)[which(colnames(cpgs) %in% names6)],"_6")
colnames(cpgs)

cpgs_overlapping=cpgs[rownames(cpgs) %in% intersecting ,]#700



####cell_composition---------------------------------
cell_com_tubingen=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction/cellcomposition.tsv',data.table=FALSE)

## without min_amx normalization

cell_com_tubingen_w=cell_com_tubingen
rownames(cell_com_tubingen_w)=phenotypes_cleaned$sample_name
rownames(cell_com_tubingen_w)[which(rownames(cell_com_tubingen_w) %in% names3)]<-paste0(rownames(cell_com_tubingen_w)[which(rownames(cell_com_tubingen_w) %in% names3)],"_3")
rownames(cell_com_tubingen_w)[which(rownames(cell_com_tubingen_w) %in% names5)]<-paste0(rownames(cell_com_tubingen_w)[which(rownames(cell_com_tubingen_w) %in% names5)],"_5")
rownames(cell_com_tubingen_w)[which(rownames(cell_com_tubingen_w) %in% names6)]<-paste0(rownames(cell_com_tubingen_w)[which(rownames(cell_com_tubingen_w) %in% names6)],"_6")
rownames(cell_com_tubingen_w)
cell_com_tubingen_w=select(cell_com_tubingen_w,-Name,-Group)
cell_com_tubingen_w_t=t(cell_com_tubingen_w)

##min_amx normalization

cell_com_tubingen$CD8T=min_max_norm(cell_com_tubingen$CD8T)
cell_com_tubingen$CD4T=min_max_norm(cell_com_tubingen$CD4T)
cell_com_tubingen$NK=min_max_norm(cell_com_tubingen$NK)
cell_com_tubingen$Bcell=min_max_norm(cell_com_tubingen$Bcell)
cell_com_tubingen$Mono=min_max_norm(cell_com_tubingen$Mono)
cell_com_tubingen$Gran=min_max_norm(cell_com_tubingen$Gran)


rownames(cell_com_tubingen)=phenotypes_cleaned$sample_name
rownames(cell_com_tubingen)[which(rownames(cell_com_tubingen) %in% names3)]<-paste0(rownames(cell_com_tubingen)[which(rownames(cell_com_tubingen) %in% names3)],"_3")
rownames(cell_com_tubingen)[which(rownames(cell_com_tubingen) %in% names5)]<-paste0(rownames(cell_com_tubingen)[which(rownames(cell_com_tubingen) %in% names5)],"_5")
rownames(cell_com_tubingen)[which(rownames(cell_com_tubingen) %in% names6)]<-paste0(rownames(cell_com_tubingen)[which(rownames(cell_com_tubingen) %in% names6)],"_6")
rownames(cell_com_tubingen)

cell_com_tubingen=select(cell_com_tubingen,-Name,-Group)
cell_com_tubingen_t=t(cell_com_tubingen)

###phenotypes

easy_measure_phenos <- phenotypes_cleaned %>%
  dplyr::select("sample_name",
                "bmi","waist","hip","crp_mg_l","glucose_000","cholesterol_mmol_l",
                "hdl_mmol_l","ldl_mmol_l","triglycerides_mmol_l","got_ast_u_l","gpt_alt_u_l","ggt_u_l","creatinine_micromol_l"
  )

pheno_easy<- easy_measure_phenos[which(easy_measure_phenos$sample_name %in% gsub("_.*","",colnames(cpgs_overlapping))),]



cluster_easy_pheno<- pheno_preper(pheno_easy,names3,names5,names6)


whole_pheno=rbind(cell_com_tubingen_t,cluster_easy_pheno)

whole_pheno_w=rbind(cell_com_tubingen_w_t,cluster_easy_pheno)

Pam_function(dataset=rbind(cpgs_overlapping,whole_pheno),name="701cpgs_Tübingen_pheno_cellcompistion")
Pam_function(dataset=rbind(cpgs_overlapping,whole_pheno_w),name="701cpgs_without_Tübingen_pheno_cellcompistion")
