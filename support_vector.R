####support vector machine
###let's try with 2020

library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidylog)
library(janitor) #clean_names(df) or %>% clean_names()
library(ragg) #use agg_png to have quicker plots
library(report) #for session report
library(tibble)
library(e1071)
library(randomForest)
library(progress)
library(ranger)
library(caret)
# options ------

#loading
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")
phenodata <- fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/raw/pheno.csv",skip=1)
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/new_combat.Rdata")

###use potsdam as test
potsdam_2020=as.data.frame(new_mycombat)


methylated_2020 <- as.data.frame(my_combat_dmp_high_frq)
rm(my_combat_dmp_high_frq)
potsdam_2020=potsdam_2020[rownames(potsdam_2020) %in% rownames(methylated_2020),]

#pheno---------------------------



phenotypes_cleaned <- phenodata %>%
  janitor::clean_names(
    case = "snake",
    replace = c("µ" = "micro", `'` = "", `"` = "", `%` = "_percent_", `#` = "_number_")
  ) %>% 
  mutate(cluster_tueftulip = case_when(
    cluster_tueftulip==3 ~ "C3",
    cluster_tueftulip==5 ~ "C5",
    cluster_tueftulip==6 ~ "C6",
  )
  ) %>%
  mutate(bmi = as.numeric(bmi)) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(waist = as.numeric(waist)) %>%
  mutate(hip = as.numeric(hip)) %>%
  mutate(crp_mg_l = as.numeric(crp_mg_l)) %>%
  mutate(hdl_mmol_l = as.numeric(hdl_mmol_l)) %>%
  mutate(ldl_mmol_l = as.numeric(ldl_mmol_l)) %>%
  mutate(cholesterol_mmol_l = as.numeric(cholesterol_mmol_l)) %>%
  mutate(got_ast_u_l = as.numeric(got_ast_u_l)) %>%
  mutate(ggt_u_l = as.numeric(ggt_u_l)) %>%
  mutate(creatinine_micromol_l= as.numeric(creatinine_micromol_l)) %>%
  mutate(mr_l_fadj = as.numeric(mr_l_fadj)) %>%
  mutate(mr_scat_l = as.numeric(mr_scat_l)) %>%
  mutate(mr_vat_l = as.numeric(mr_vat_l)) %>%
  mutate(di_calc = as.numeric(di_calc)) %>%
  mutate(isi_calc = as.numeric(isi_calc)) %>%
  dplyr::select(-sample_group,-pool_id,-sample_plate,-sample_well,-sentrix_id,-sentrix_position,-prs)

rownames(phenotypes_cleaned)<-phenotypes_cleaned$sample_name

names3 <- as.character(phenotypes_cleaned$sample_name[which(phenotypes_cleaned$cluster_tueftulip == "C3")])
names5 <- as.character(phenotypes_cleaned$sample_name[which(phenotypes_cleaned$cluster_tueftulip == "C5")])
names6 <- as.character(phenotypes_cleaned$sample_name[which(phenotypes_cleaned$cluster_tueftulip == "C6")])

#
colnames(methylated_2020)=gsub('_3','',colnames(methylated_2020))
colnames(methylated_2020)=gsub('_5','',colnames(methylated_2020))
colnames(methylated_2020)=gsub('_6','',colnames(methylated_2020))


methylated_2020_2<-copy(methylated_2020)

##i don't want to split data
####just to add coulmn cluster in my data
methylated_2020_2<-data.frame(t(methylated_2020))
methylated_2020_2$group<-"NA"
methylated_2020_2$group[rownames(methylated_2020_2) %in% names3] <- "C3"
methylated_2020_2$group[rownames(methylated_2020_2) %in% names5] <- "C5"
methylated_2020_2$group[rownames(methylated_2020_2) %in% names6] <- "C6"

###dont forget it mustr be factor 
methylated_2020_2$group=as.factor(methylated_2020_2$group)


pheno_pot=fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno.csv",skip=1,data.table=FALSE)

names3_pot<- pheno_pot$Sample_Name[which(pheno_pot$cluster_tueftulip == 3)]
names5_pot<- pheno_pot$Sample_Name[which(pheno_pot$cluster_tueftulip == 5)]
names6_pot<- pheno_pot$Sample_Name[which(pheno_pot$cluster_tueftulip == 6)]


potsdam_2020_2<-data.frame(t(potsdam_2020))
potsdam_2020_2$group<-"NA"
potsdam_2020_2$group[rownames(potsdam_2020_2) %in% names3_pot] <- "C3"
potsdam_2020_2$group[rownames(potsdam_2020_2) %in% names5_pot] <- "C5"
potsdam_2020_2$group[rownames(potsdam_2020_2) %in% names6_pot] <- "C6"

potsdam_2020_2$group=as.factor(potsdam_2020_2$group)

##change labeling
methylation_2=methylated_2020[rownames(methylated_2020) %in% rownames(potsdam_2020),]
methylation_22<-data.frame(t(methylation_2))
methylation_22$group<-"NA"
methylation_22$group[rownames(methylation_22) %in% names3] <- "C3"
methylation_22$group[rownames(methylation_22) %in% names5] <- "C5"
methylation_22$group[rownames(methylation_22) %in% names6] <- "C6"


###dont forget it must be factor 
methylation_22$group=as.factor(methylation_22$group)

###train control
trctrl <- trainControl(method = "LOOCV", number = 20, repeats = 3,verboseIter = T)

methyl_tub_train=train(group~., data = methylated_2020_2, method = 'svmPoly',trControl=trctrl,tuneLength = 10)




methyl_pot_train=train(group~., data = potsdam_2020_2,method='svmPoly',trControl=trctrl,preProcess = c("center", "scale"),tuneLength = 10)


save(methyl_tub_train,file='support_tubingen.Rdata')
save(methyl_pot_train,file='support_potsdam.Rdata')
