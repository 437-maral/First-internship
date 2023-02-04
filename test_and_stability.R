
####################################################################################################

#####################################################################################################

# libraries ----------------------------------------------------------------
library(cluster)
library(factoextra)
library(M3C)
library(openxlsx)
library(cvms)
library(rgl)
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(reshape2)
library(ggplot2)
library(progress)


# setwd -------------------------------------------------------------------

setwd('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/two_cohorts')


# readin ------------------------------------------------------------------

pheno=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/two_cohorts/raw/Thuebingen_Potsdam_cohorts_2outlier_removed.csv',skip = 1,data.table = FALSE)
load('dmp_twocohort_potsdam_2outliersrm.RData')
load('two_combats_2outliersrm.Rdata')


# prep pheno --------------------------------------------------------------


pheno$Sample_Name <- as.character(pheno$Sample_Name)#like last time 
pheno$cluster_tueftulip <- as.factor(pheno$cluster_tueftulip)

pheno_potsdam<-pheno %>%
  filter(Sample_Plate=="WG7035812-MSA4")



# prep data ---------------------------------------------------------------

######first start with dmps to find a subset size 
####### I chose  dmps in intersecting region , because the have high acuuracy 


dmp6_5=dmptest_potsdam$`6_to_5`[dmptest_potsdam$`6_to_5`$P.Value<0.05,] #40417
dmp3_6=dmptest_potsdam$`6_to_3`[dmptest_potsdam$`6_to_3`$P.Value<0.05,] #42914
dmp5_3=dmptest_potsdam$`5_to_3`[dmptest_potsdam$`5_to_3`$P.Value<0.05,] #22784

dmp6_5=dmptest_potsdam$`6_to_5`[dmptest_potsdam$`6_to_5`$P.Value<0.01,] # 8073
dmp3_6=dmptest_potsdam$`6_to_3`[dmptest_potsdam$`6_to_3`$P.Value<0.01,] #3540
dmp5_3=dmptest_potsdam$`5_to_3`[dmptest_potsdam$`5_to_3`$P.Value<0.01,] #9788

dmp6_5_adj=dmp6_5[dmp6_5$adj.P.Val<0.05,] #nothing
dmp3_6_adj=dmp6_5[dmp6_5$adj.P.Val<0.05,] #nothing
dmp5_3_adj=dmp6_5[dmp6_5$adj.P.Val<0.05,] #nothing

##add Name column for cpgs
dmp6_5$Name=rownames(dmp6_5)
dmp3_6$Name=rownames(dmp3_6)
dmp5_3$Name=rownames(dmp5_3)


dmps <- c(dmp6_5$Name,dmp5_3$Name,dmp3_6$Name)




Potsdam_cohort=two_myCombats[,colnames(two_myCombats) %in% pheno_potsdam$Sample_Name]

######all_dmps

mycombat_dmps=Potsdam_cohort[rownames(Potsdam_cohort) %in% dmps ,]#   46  /19020
#now we realy only have each cpg one time!


mycombat_dmps=data.frame(mycombat_dmps)


#needs to be done, otherwise our sample names would be taken as integer value
names3 <- paste0("X",pheno_potsdam$Sample_Name[which(pheno_potsdam$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno_potsdam$Sample_Name[which(pheno_potsdam$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno_potsdam$Sample_Name[which(pheno_potsdam$cluster_tueftulip == 6)])

colnames(mycombat_dmps)[which(colnames(mycombat_dmps) %in% names3)]<-paste0(colnames(mycombat_dmps)[which(colnames(mycombat_dmps) %in% names3)],"_3")
colnames(mycombat_dmps)[which(colnames(mycombat_dmps) %in% names5)]<-paste0(colnames(mycombat_dmps)[which(colnames(mycombat_dmps) %in% names5)],"_5")
colnames(mycombat_dmps)[which(colnames(mycombat_dmps) %in% names6)]<-paste0(colnames(mycombat_dmps)[which(colnames(mycombat_dmps) %in% names6)],"_6")
colnames(mycombat_dmps)


####intersecting_dmps

dmps <- c(dmp6_5$Name,dmp5_3$Name,dmp3_6$Name)
duplicated  <- dmps[duplicated(dmps)] #22k



mycombat_intersecting=Potsdam_cohort[rownames(Potsdam_cohort) %in% duplicated ,]#2379   46

mycombat_intersecting=data.frame(mycombat_intersecting)


names3_int<- paste0("X",pheno_potsdam$Sample_Name[which(pheno_potsdam$cluster_tueftulip == 3)])
names5_int<- paste0("X",pheno_potsdam$Sample_Name[which(pheno_potsdam$cluster_tueftulip == 5)])
names6_int<- paste0("X",pheno_potsdam$Sample_Name[which(pheno_potsdam$cluster_tueftulip == 6)])

colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names3_int)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names3_int)],"_3")
colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names5_int)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names5_int)],"_5")
colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names6_int)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names6_int)],"_6")
colnames(mycombat_intersecting)







# get sample size with max acc --------------------------------------------

####here ,we chose 190 times , becuase we had 19020 cpgs
experiment_results <- list()
size <- 1
size_vector <- c()
acc <- c()
pb <-progress_bar$new(total=190,clear=FALSE,format= "(:spin) [:bar] :current/:total (:percent) in :elapsedfull")
for (i in 1:190){ 
  sample <- mycombat_dmps[sample(nrow(mycombat_dmps), size=size,replace=FALSE), ]
  pb$tick()
  pa <-pam(t(sample), 3, diss =  FALSE,
           metric = "manhattan", 
           medoids =  "random",
           nstart =  10 ,
           stand = FALSE, cluster.only = FALSE,
           do.swap = TRUE,
           keep.diss = FALSE,
           keep.data = FALSE,
           variant = "faster",
           trace.lev = 0)
  
  pa$data <- t(sample)
  clustering <- as.data.frame(pa[["clustering"]])
  clustering <- rownames_to_column(clustering)
  #get only the number of the cluster
  clustering$rowname <- gsub(".+_","", clustering$rowname) %>% as.factor()
  colnames(clustering) <-c("target","prediciton")
  #important to change 3 first 
  #sofyas clustering
  clustering_h<-copy(clustering)
  clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
  clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
  d <- table(clustering_h)
  sofya <- sum(diag(d))/sum(d) #overall accuracy
  
  clustering_h<-copy(clustering)
  clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
  d <- table(clustering_h)
  markus <- sum(diag(d))/sum(d) #overall accuracy
  
  clustering_h<-copy(clustering)
  clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
  clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton)
  d <- table(clustering_h)
  uk1 <- sum(diag(d))/sum(d) #overall accuracy
  
  clustering_h<-copy(clustering)
  clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
  clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton) 
  d <- table(clustering_h)
  uk2 <- sum(diag(d))/sum(d) #overall accuracy
  
  clustering_h<-copy(clustering)
  clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton)
  clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton)  
  d <- table(clustering_h)
  uk3 <- sum(diag(d))/sum(d) #overall accuracy
  
  clustering_h<-copy(clustering)
  clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
  clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton)
  clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
  d <- table(clustering_h)
  uk4 <- sum(diag(d))/sum(d) #overall accuracy
  
  acc[i] <- max(sofya,markus,uk1,uk2,uk3,uk4)
  size_vector[i] <- size
  size <- size +100
}



save(acc,file='acc_test.Rdata')
load('acc_test.Rdata')

df <- as.data.frame(cbind(acc,size_vector))

plot1 <- ggplot(df, aes(size_vector)) +
  geom_line(aes(y = acc)) +
  labs(x="subset Size")+ 
  ggtitle("Accuracy of PAM clustering with Top n DMPs") +
  theme_grey(base_size = 15)+
  geom_hline(yintercept =0.80,color='red')
plot1

png("first_check.png")
plot1
dev.off()


max(df$acc) #0.869
df[which(df$acc ==max(df$acc)),] #1801 subset size

fwrite(df,"dataframe_for_image.csv")

maximum_acc=fread("dataframe_for_image.csv")
# check the stability -----------------------------------------------------

gc()

#check the stability
experiment_results <- list()
pb <-progress_bar$new(total=100,clear=FALSE,format= "(:spin) [:bar] :current/:total (:percent) in :elapsedfull")

for (j in 1:100){
  pb$tick()
  
  size <- 1
  size_vector <- c()
  acc <- c()
  for (i in 1:100){
    
    sample <- mycombat_dmps[sample(nrow(mycombat_dmps), size=size), ] 
    pa <-pam(t(sample), 3, diss =  FALSE,
             metric = "manhattan", 
             medoids =  "random",
             nstart =  10 ,
             stand = FALSE, cluster.only = FALSE,
             do.swap = TRUE,
             keep.diss = FALSE,
             keep.data = FALSE,
             variant = "faster",
             trace.lev = 0)
    
    pa$data <- t(sample)
    clustering <- as.data.frame(pa[["clustering"]])
    clustering <- rownames_to_column(clustering)
    #get only the number of the cluster
    clustering$rowname <- gsub(".+_","", clustering$rowname) %>% as.factor()
    colnames(clustering) <-c("target","prediciton")
    #important to change 3 first 
    #sofyas clustering
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
    d <- table(clustering_h)
    sofya <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
    d <- table(clustering_h)
    markus <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton)
    d <- table(clustering_h)
    uk1 <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton) 
    d <- table(clustering_h)
    uk2 <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton)  
    d <- table(clustering_h)
    uk3 <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
    d <- table(clustering_h)
    uk4 <- sum(diag(d))/sum(d) #overall accuracy
    
    acc[i] <- max(sofya,markus,uk1,uk2,uk3,uk4)
    size_vector[i] <- size
    size <- size +100
  }
  experiment_results[[j]] <-acc 
  
}

save(experiment_results,file='acc_stability.Rdata')


experiment_results_as_mat <- matrix(unlist(experiment_results), ncol=100,byrow = T)
colnames(experiment_results_as_mat) <- size_vector #cpgs in sample
rownames(experiment_results_as_mat) <- c(1:100) #experiemnt nr

#plot experiment (10 times repeat) results
#melt data 
experiment_results_as_mat_melted <-melt(experiment_results_as_mat)
#smoothed plot
My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18),
  plot.title = element_text(size = 22))

plot2 <- ggplot(experiment_results_as_mat_melted , aes(Var2, value)) +  
  geom_point() +
  geom_smooth(se = TRUE,method = "loess")+
  xlab("Subset size") + ylab("Accuracy")+
  ggtitle("Accuracy of PAM with random subsets of sizes ")+
  theme_light()
png("stability.png")
plot2 + My_Theme
dev.off()

fwrite(experiment_results_as_mat,"dataframe_of_stability.csv")






# do the wrapper ----------------------------------------------------------


max_acc_total<-0
experiment_results <- list()
acc <- list()
for (j in 1:10){
  cat(paste0("starting run:",j,"\n"))
  #loop over pam acc for a sample subset of 1000 CpG and save those subsets that reach an Acc over 85%
  acc[[j]] <- list()
  subsets_good_acc_85 <- list()
  pb <-progress_bar$new(total=100000,clear=FALSE,
                        format= "(:spin) :what [:bar] :current/:total (:percent) in :elapsedfull")
  for (i in 1:100000){
    acc[[j]][[i]]<-list()
    pb$tick(tokens = list(what= j))
    part <- mycombat_dmps[sample(nrow(mycombat_dmps), size=2000), ] 
    pa <-pam(t(part), 3, diss =  FALSE,
             metric = "manhattan", 
             medoids =  "random",
             nstart =  10 ,
             stand = FALSE, cluster.only = FALSE,
             do.swap = TRUE,
             keep.diss = FALSE,
             keep.data = FALSE,
             variant = "faster",
             trace.lev = 0)
    
    pa$data <- t(part)
    clustering <- as.data.frame(pa[["clustering"]])
    clustering <- rownames_to_column(clustering)
    #get only the number of the cluster
    clustering$rowname <- gsub(".+_","", clustering$rowname) %>% as.factor()
    colnames(clustering) <-c("target","prediciton")
    #important to change 3 first because it is in both target and prediciton
    
    #####markus put this in
    
    #sofyas clustering
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["sofya"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["markus"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["uk1"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton) 
    d <- table(clustering_h)
    acc[[j]][[i]][["uk2"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton)  
    d <- table(clustering_h)
    acc[[j]][[i]][["uk3"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["uk4"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    max_acc <-max(acc[[j]][[i]][["sofya"]], acc[[j]][[i]][["markus"]],
                  acc[[j]][[i]][["uk1"]],acc[[j]][[i]][["uk2"]],
                  acc[[j]][[i]][["uk3"]],acc[[j]][[i]][["uk4"]])
    
    if(max_acc > max_acc_total){
      cat(paste0("new max: ", max_acc,"\n"))
      max_acc_total<-max_acc
    }
    
    if (max_acc > 0.85){
      subsets_good_acc_85[[i]] <- rownames(part)
      names(subsets_good_acc_85[[i]])<-names(acc[[j]][[i]])[which.min(acc[[j]][[i]])]
    }

  }
  experiment_results[[j]] <-subsets_good_acc_85 
}


save(experiment_results, file= "10times_Potsdam_2000Samples_subsets.RData")
save(acc, file = "10times_Potsdam_2000Samples_subsets_acc.RData")



 #### work on tÜbingen cohort


### conider all dmps 



load('dmp_twocohort_tubingen_2outliersrm.Rdata')

dmp6_5_tue=dmptest_tubingen$`6_to_5`[dmptest_tubingen$`6_to_5`$P.Value<0.05,]
dmp3_6_tue=dmptest_tubingen$`3_to_6`[dmptest_tubingen$`3_to_6`$P.Value<0.05,]
dmp3_5_tue=dmptest_tubingen$`3_to_5`[dmptest_tubingen$`3_to_5`$P.Value<0.05,]


##add Name column for cpgs
dmp6_5_tue$Name=rownames(dmp6_5_tue)
dmp3_6_tue$Name=rownames(dmp3_6_tue)
dmp3_5_tue$Name=rownames(dmp3_5_tue)

library(tidyverse)
pheno_Tubingen=pheno %>%
  filter (str_detect(Sample_Plate,'A4726'))

#####pam

Tubingen_cohort=two_myCombats[,colnames(two_myCombats) %in% pheno_Tubingen$Sample_Name]

all_cpgs_tu=c(dmp6_5_tue$Name,dmp3_6_tue$Name,dmp3_5_tue$Name)#65895

mycombat_alldmps_tue=Tubingen_cohort[rownames(Tubingen_cohort) %in% all_cpgs_tu,]#65895    72
mycombat_alldmps_tue=data.frame(mycombat_alldmps_tue)

names3_tue<- paste0("X",pheno_Tubingen$Sample_Name[which(pheno_Tubingen$cluster_tueftulip == 3)])
names5_tue<- paste0("X",pheno_Tubingen$Sample_Name[which(pheno_Tubingen$cluster_tueftulip == 5)])
names6_tue<- paste0("X",pheno_Tubingen$Sample_Name[which(pheno_Tubingen$cluster_tueftulip == 6)])


colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names3_tue)]<-paste0(colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names3_tue)],"_3")
colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names5_tue)]<-paste0(colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names5_tue)],"_5")
colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names6_tue)]<-paste0(colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names6_tue)],"_6")
colnames(mycombat_alldmps_tue)


max_acc_total<-0
experiment_results <- list()
acc <- list()
for (j in 1:10){
  cat(paste0("starting run:",j,"\n"))
  #loop over pam acc for a sample subset of 1000 CpG and save those subsets that reach an Acc over 85%
  acc[[j]] <- list()
  subsets_good_acc_85 <- list()
  pb <-progress_bar$new(total=100000,clear=FALSE,
                        format= "(:spin) :what [:bar] :current/:total (:percent) in :elapsedfull")
  for (i in 1:100000){
    acc[[j]][[i]]<-list()
    pb$tick(tokens = list(what= j))
    part <- mycombat_alldmps_tue[sample(nrow(mycombat_alldmps_tue), size=1000), ] 
    pa <-pam(t(part), 3, diss =  FALSE,
             metric = "manhattan", 
             medoids =  "random",
             nstart =  10 ,
             stand = FALSE, cluster.only = FALSE,
             do.swap = TRUE,
             keep.diss = FALSE,
             keep.data = FALSE,
             variant = "faster",
             trace.lev = 0)
    
    pa$data <- t(part)
    clustering <- as.data.frame(pa[["clustering"]])
    clustering <- rownames_to_column(clustering)
    #get only the number of the cluster
    clustering$rowname <- gsub(".+_","", clustering$rowname) %>% as.factor()
    colnames(clustering) <-c("target","prediciton")
    #important to change 3 first because it is in both target and prediciton
    
    #####markus put this in
    
    #sofyas clustering
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["sofya"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(1,3,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["markus"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,6,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["uk1"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton) 
    d <- table(clustering_h)
    acc[[j]][[i]][["uk2"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,5,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(3,6,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(2,3,clustering_h$prediciton)  
    d <- table(clustering_h)
    acc[[j]][[i]][["uk3"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    clustering_h<-copy(clustering)
    clustering_h$prediciton <- gsub(1,6,clustering_h$prediciton) 
    clustering_h$prediciton <- gsub(2,5,clustering_h$prediciton)
    clustering_h$prediciton <- gsub(3,3,clustering_h$prediciton)
    d <- table(clustering_h)
    acc[[j]][[i]][["uk4"]] <- sum(diag(d))/sum(d) #overall accuracy
    
    max_acc <-max(acc[[j]][[i]][["sofya"]], acc[[j]][[i]][["markus"]],
                  acc[[j]][[i]][["uk1"]],acc[[j]][[i]][["uk2"]],
                  acc[[j]][[i]][["uk3"]],acc[[j]][[i]][["uk4"]])
    
    if(max_acc > max_acc_total){
      cat(paste0("new max: ", max_acc,"\n"))
      max_acc_total<-max_acc
    }
    
    if (max_acc > 0.80){
      subsets_good_acc_85[[i]] <- rownames(part)
      names(subsets_good_acc_85[[i]])<-names(acc[[j]][[i]])[which.min(acc[[j]][[i]])]
    }
    
  }
  experiment_results[[j]] <-subsets_good_acc_85 
}





load('markus_dec2022/10times_tubingen_1000Samples_subsets.RData')
load('markus_dec2022/10times_tubingen_1000Samples_subsets_acc.RData')

df=readRDS('markus_dec2022/tabled_freq_sorted.rds')















#####just to check whether this has high accuracy or not
load("markus_dec2022/10times_Potsdam_2000Samples_subsets.RData")



subsets_clearedf_85 <- unlist(experiment_results[!unlist(lapply(experiment_results,is.null))])
tabled_freq <- as.data.frame(table(subsets_clearedf_85))
tabled_freq_sorted <- tabled_freq[order(-tabled_freq$Freq), ]
cpgs_high_freq_acc <- tabled_freq %>%  filter(Freq>600) %>% pull(subsets_clearedf_85)


###I want to compare with whole cohort
###let's overalpp with intersecting dmp 

my_combat_dmp_high_frq_int<- subset(mycombat_intersecting,rownames(mycombat_intersecting) %in% cpgs_high_freq_acc)


my_combat_dmp_high_frq <- subset(Potsdam_cohort,rownames(Potsdam_cohort) %in% cpgs_high_freq_acc) #2819

# do pam 




