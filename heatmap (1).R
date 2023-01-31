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
###
source("L:/!DIAB/Masoumeh Vojood/functions/min_max_norm.R")
source("L:/!DIAB/Masoumeh Vojood/functions/PAM.R")
source('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/codes/pam_euclidean.R')
######directory 
setwd("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/")
load('Combat.Rdata')

load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")
pheno<-fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno1.csv",skip=1,data.table=FALSE)


####before that you should run this script dmp_overlapping_with_2020_check.R 


cpgs_new1 <- calcuform_meth(as.data.frame(cpgs2020_mean_helper[,1:4]))
cpgs_new1 <- na.omit(cpgs_new1)



cpgs_new1= data.frame(cpgs_new1 )%>% mutate(group=cpgs2020_mean_helper$group)
#####
#filtering

#only 3_5
cpgs_new1=data.frame(cpgs_new1)

##first comparison 3,5
high_35= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 >= 1 & AVG_5  <= -0.4 & AVG_6  <= -0.3)
high_35=high_35[high_35$group=="only_3_5",]
dim(high_35) 


high_53= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <=-0.4 & AVG_5  >= 1 & AVG_6  <= -0.3)
high_53=high_53[high_53$group=="only_3_5",]
dim(high_53) 


#### first comparison 3,6

high6_3= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 >= 1 & AVG_5  <= -0.3 & AVG_6  <= -0.4)
high6_3=high6_3[high6_3$group=="only_3_6",]
dim(high6_3) 



high6_6= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <=-0.4 & AVG_5  <= -0.3 & AVG_6  >=1)
high6_6=high6_6[high6_6$group=="only_3_6",]
dim(high6_6) 



###### for comparison 5,6

high5_6= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <=-0.3 & AVG_5  >=1 & AVG_6 <= -0.4)
high5_6=high5_6[high5_6$group=="only_5_6",]
dim(high5_6) 


high5_6_6= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <=-0.3 & AVG_5 <= -0.4 & AVG_6 >=1)
high5_6_6=high5_6_6[high5_6_6$group=="only_5_6",]
dim(high5_6_6) 


#####for comparison 35_36
high_3_35_36= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 >=1,AVG_5 <= -0.4 & AVG_6  <= -0.4)
high_3_35_36=high_3_35_36[high_3_35_36$group=='only_35_36',]
dim(high_3_35_36) 

####for comparison 35_36
low_3_35_36= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <= -1 & AVG_5 >=0.5 & AVG_6 >= 0.5)
low_3_35_36=low_3_35_36[low_3_35_36$group=="only_35_36",]
dim(low_3_35_36)




#### for comparison 35_56
high_5_35_56= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <= -0.4 & AVG_5 >1 & AVG_6 <= -0.4)
high_5_35_56=high_5_35_56[high_5_35_56$group=="only_35_56",]
dim(high_5_35_56)

low_5_35_56= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 >= 0.5 & AVG_5 <= -1 & AVG_6 >= 0.5)
low_5_35_56=low_5_35_56[low_5_35_56$group=="only_35_56",]
dim(low_5_35_56)

###for comparison 36_56

high_6_36_56= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 <= -0.4 & AVG_5 <= -0.4  & AVG_6 >= 1)
high_6_36_56=high_6_36_56[high_6_36_56$group=="only_36_56",]
dim(high_6_36_56)


low_6_36_56= cpgs_new1 %>% select(AVG_3,group,AVG_5, AVG_6) %>% filter(AVG_3 >=0.5 & AVG_5 >=0.5  & AVG_6  <= -1)
low_6_36_56=low_6_36_56[low_6_36_56$group=="only_36_56",]
dim(low_6_36_56)




new1=rbind(high_35,high_53,high6_3,high6_6,high5_6,high5_6_6,high_3_35_36,low_3_35_36,low_5_35_56,high_6_36_56,low_6_36_56,high_5_35_56)

all_cpgs=rownames(new1)

cpgs_disc<- my_combat_dmp_high_frq%>% filter(rownames(.) %in% all_cpgs) #450

pam_function(input=cpgs_disc,class="sofya") #81.94
load('Combat.Rdata')
com=data.frame(myCombat)

rep2020<-com %>% filter(rownames(.) %in% all_cpgs) #450



names3 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 6)])


colnames(rep2020)[which(colnames(rep2020) %in% names3)]<-paste0(colnames(rep2020)[which(colnames(rep2020) %in% names3)],"_3")
colnames(rep2020)[which(colnames(rep2020) %in% names5)]<-paste0(colnames(rep2020)[which(colnames(rep2020) %in% names5)],"_5")
colnames(rep2020)[which(colnames(rep2020) %in% names6)]<-paste0(colnames(rep2020)[which(colnames(rep2020) %in% names6)],"_6")
colnames(rep2020)

pam_function(input=rep2020,class="unknown4")#41.67

pam__euclidean(input=rep2020,class="unknown4")#41.676
###try
pam_function(input=rep2020,class="sofya") #29
pam_function(input=rep2020,class="markus")#32
pam_function(input=rep2020,class="unknown1")#31
pam_function(input=rep2020,class="unknown2")#33
pam_function(input=rep2020,class="unknown3") #29
pam_function(input=rep2020,class="unknown4")#41


######heatmap for 450 cpgs

com=data.frame(myCombat)
rep450<-com %>% filter(rownames(.) %in% all_cpgs) #450
rep450$ID=rownames(rep450)

######I want to show how 450 cpgs show us unique methylation



calcuform_meth <- function(df){
  rownames(df)<-make.names(df[,1],unique=TRUE)
  df[,1]<-NULL
  log.z <- t(scale(t(df)))  #don't do the log2 calculation
  
  mat = as.matrix(log.z) #matrix
  rownames(mat) = rownames(log.z) #rownames
  colnames(mat) = colnames(log.z) #colnames
  
  return(mat)
}


####colnames are just samples , so we need to rename them regarding to cluster's name
names3 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno$Sample_Name[which(pheno$cluster_tueftulip == 6)])


colnames(rep450)[which(colnames(rep450) %in% names3)]<-paste0(colnames(rep450)[which(colnames(rep450) %in% names3)],"_3")
colnames(rep450)[which(colnames(rep450) %in% names5)]<-paste0(colnames(rep450)[which(colnames(rep450) %in% names5)],"_5")
colnames(rep450)[which(colnames(rep450) %in% names6)]<-paste0(colnames(rep450)[which(colnames(rep450) %in% names6)],"_6")
colnames(rep450)



######getting mean of each cluster



df33<- rep450 %>%
  dplyr::select(matches("_3$")) %>%
  rowwise() %>%
  mutate(AVG_3= mean(c_across())) 
df55<- rep450 %>%
  dplyr::select(matches("_5$")) %>%
  rowwise() %>%
  mutate(AVG_5= mean(c_across()))
df66<- rep450 %>%
  dplyr::select(matches("_6$")) %>%
  rowwise() %>%
  mutate(AVG_6= mean(c_across()))                  


cpg450<- data.frame(cbind(ID=rep450$ID,df33,df55,df66)) %>%
  select(ID,AVG_3,AVG_5,AVG_6)


preped450 <- calcuform_meth(as.data.frame(cpg450))
preped450 <- na.omit(preped450)
col_fun450= colorRamp2(seq(min(preped450[is.finite(preped450)]), max(preped450[is.finite(preped450)]), length = 3), c("#0088DD", "white", "#EE0000"))





Heatmap(
  name="mean methylation",
  matrix=preped450,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_row_dend=FALSE,
  row_dend_width = unit(4, "cm"),
  col = col_fun450,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  row_title_rot = 0,
  #column_labels= c(rep("",ncol(prepared_sorted))),
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 3),
  #border_gp = gpar(col = "#636363", lwd = 1,lty=2),
  column_title=paste0( nrow(preped450)," CpGs "),
  column_title_side="bottom"
)








