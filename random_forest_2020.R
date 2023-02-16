library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidylog)
library(janitor) #clean_names(df) or %>% clean_names()
library(ragg) #use agg_png to have quicker plots
library(report) #for session report
library(tibble)
library(randomForest)
library(progress)
library(ranger)
library(e1071)
# options -----------------------------------------------------------------
"%!in%" <- function(x, y) !("%in%"(x, y))



# readin ------------------------------------------------------------------

load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")
phenodata <- fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/raw/pheno.csv",skip=1)
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/new_combat.Rdata")

###use potsdam as test
potsdam_2020=as.data.frame(new_mycombat)


methylated_2020 <- as.data.frame(my_combat_dmp_high_frq)
rm(my_combat_dmp_high_frq)
potsdam_2020=potsdam_2020[rownames(potsdam_2020) %in% rownames(methylated_2020),]

####




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

###
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


#potsdam manipulation 
###look at  potsdam phenotypes
pheno_pot=fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno.csv",skip=1,data.table=FALSE)

names3_pot<- pheno_pot$Sample_Name[which(pheno_pot$cluster_tueftulip == 3)]
names5_pot<- pheno_pot$Sample_Name[which(pheno_pot$cluster_tueftulip == 5)]
names6_pot<- pheno_pot$Sample_Name[which(pheno_pot$cluster_tueftulip == 6)]


potsdam_2020_2<-data.frame(t(potsdam_2020))
potsdam_2020_2$group<-"NA"
potsdam_2020_2$group[rownames(potsdam_2020_2) %in% names3_pot] <- "C3"
potsdam_2020_2$group[rownames(potsdam_2020_2) %in% names5_pot] <- "C5"
potsdam_2020_2$group[rownames(potsdam_2020_2) %in% names6_pot] <- "C6"

###dont forget it mustr be factor 
potsdam_2020_2$group=as.factor(potsdam_2020_2$group)

###change the labeling
###becuase we have different cpg for both data set, i used potsdam cpgs as reference for tübingen 
methylation_2=methylated_2020[rownames(methylated_2020) %in% rownames(potsdam_2020),]
methylation_22<-data.frame(t(methylation_2))
methylation_22$group<-"NA"
methylation_22$group[rownames(methylation_22) %in% names3] <- "C3"
methylation_22$group[rownames(methylation_22) %in% names5] <- "C5"
methylation_22$group[rownames(methylation_22) %in% names6] <- "C6"

###dont forget it mustr be factor 
methylation_22$group=as.factor(methylation_22$group)

###just to look at cofusion matrix

model_22<- ranger(group ~ ., 
                data = methylation_22, 
                num.trees = 1000,
                save.memory=TRUE,
                importance = "impurity_corrected") 


predTest <- predict(model_22,data= potsdam_2020_2)

table(potsdam_2020_2$group, predTest$predictions)


####to check their confusion matrix

model <- ranger(group ~ ., 
                data = methylation_22, 
                num.trees = 1000,
                save.memory=TRUE,
                importance = "impurity_corrected") 


model_1<- ranger(group ~ ., 
                data = potsdam_2020_2, 
                num.trees = 1000,
                save.memory=TRUE,
                importance = "impurity_corrected") 





#------------------lasso
best_a<-0
lowest_a<-100
n_runs <- 100


model=list()
predTest=list()

for (j in 1:n_runs) {
  pb <- progress_bar$new(
    total = n_runs, clear = FALSE, show_after = 0,
    format = "(run# :what  best# :best worst# :lowest) [:bar] :current/:total (:percent) in :elapsedfull"
  )
  pb$tick(0)
  model[[j]]<- ranger(group ~ ., 
                  data = methylation_22, 
                  num.trees = 1000,
                  save.memory=TRUE,
                  importance = "impurity_corrected") 
  # check if the model performs good enough
  
  predTest[[j]]<- predict(model[[j]], data=potsdam_2020_2)
  # save the best A
  if (mean(predTest[[j]]$predictions == potsdam_2020_2$group) > best_a) {
    best_a <- mean(predTest[[j]]$predictions == potsdam_2020_2$group)
  }
  # save the lowest A
  if (mean(predTest[[j]]$predictions == potsdam_2020_2$group) < lowest_a) {
    lowest_a <- mean(predTest[[j]]$predictions== potsdam_2020_2$group)
  }
  
  # save only if model is good
  if (mean(predTest[[j]]$predictions == potsdam_2020_2$group) > 0.7) { # ~70% have to be correct
    f<- data.frame(importance_pvalues(model[[j]], method = "janitza"))
    f$cg <- rownames(f)
    f <- f %>%
      dplyr::filter(pvalue < 0.05) 
  }
  pb$tick(tokens = list(run= j,best=round(best_a,2),lowest=round(lowest_a,2)))
}



# test --------------------------------------------------------------------
# tic()
# modelrandomforest <- randomForest(group ~ ., data = setlist_train[[1]], ntree = 100,importance = TRUE)  #68 sec
# toc()
# tic()
# modelraner <- ranger(group ~ ., data = setlist_train[[1]], num.trees = 1000) #3 sec
# toc()
# tic()
# modelrangerlowmem <- ranger(group ~ ., data = setlist_train[[1]], num.trees = 1000,save.memory=TRUE) #3sec
# toc()

ranger_full_train <-filtered_methylation[,which(colnames(filtered_methylation) %in% train)]
ranger_full_train<-data.frame(t(ranger_full_train))
#no chance..would need more then 400GB of ram!

tic()
rrf <- ranger(group ~ ., 
              data = setlist_train[[1]], 
              num.trees = 1000,
              save.memory=TRUE,
              importance = "impurity_corrected") 
toc()
rrf.prediction = predict(rrf, data=setlist_test[[1]])


rrf$confusion.matrix

impo<-data.frame(importance_pvalues(rrf, method = "janitza"))
#filter for pvalue < 0.05 

rrf.prediction = predict(rrf, data=setlist_test[[1]])
rrf.prediction

table(setlist_test[[1]]$group, predictions(rrf.prediction))




# session report ----------------------------------------------------------
writeLines(c(report_system(),"\n",report_packages(),"\n",report(sessionInfo())),sep = "\n","session_report.txt")


####################################################################