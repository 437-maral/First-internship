########I started with tubingen Cohort
library(archetypes)
library(data.table)


####set directory



load("two_combats_2outliersrm.Rdata")
load('2020cpgs.Rdata')
source('min_max_norm.R')

pheno1=fread('pheno.csv',skip = 1,data.table = FALSE)

obs_cohort=two_myCombats[, colnames(two_myCombats) %in% pheno1$Sample_Name]
tübingen_cohort=obs_cohort[rownames(obs_cohort) %in% rownames(my_combat_dmp_high_frq),]
tübingen_cohort_t=t(tübingen_cohort)

####I selected robust archetypes ,because it is more robust to outlier
Tubingen_arch<-stepArchetypes(t(tübingen_cohort), k=1:10,nrep = 100,family = archetypesFamily("robust",zalphasfn = archetypes:::ginv.zalphasfn))
save(Tubingen_arch,file='archetypes_tub')
screeplot(Tubingen_arch,type='barplot')

#####stability
# with different K
submodel=list()
for (k in 1:10) {
  for (i in 1:100) {
    subsample <- subsample <- sample(1:nrow(tübingen_cohort_t), size = 0.9 * nrow(tübingen_cohort_t))
    submodel[[k]] <- archetypes(tübingen_cohort_t[subsample, ], k = k,family = archetypesFamily("robust",zalphasfn = archetypes:::ginv.zalphasfn))
    }
}
##save model
save(submodel,file='submodel_archetypes.Rdata')


###refrence cluster
#### I am not sure which cluster refrence I should use ??????*******

### so go for 3,5,6 as cluster reference 
label_original=select(pheno1,cluster_tueftulip)
label_original=label_original[subsample,]
label_original=as.vector(label_original)

#loop
cutoffs <- seq(0, 1, 0.05)
sample_data=list()
rand=list()
for(j in 1:length(cutoffs)){
  sample_data[[j]]=list()
  rand[[j]]=list()
  for(m in 1:10){
  ###change a soft clustering to a hard clustering for rand index
  sample_data[[j]][[m]]=list()
  rand[[j]][[m]]=list()
  #submodel[[m]]$alphas for membershup cutoff
  sample_data[[j]][[m]]=apply(submodel[[m]]$alphas, 1, function(x) which.max(x >= cutoffs[j] & x < cutoffs[j+1]))
  #I got NaN for adj.rand.index functiom , so I chose rand.index function !
  rand[[j]][[m]]=rand.index(sample_data[[j]][[m]],label_original)
  }
}
####do plot





#-----------------------PCA

#####I'll go with four archetypes!!!!! here , I am not sure about the number of archetypes 
#library
library(dplyr)
library(factoextra)
library(tidyverse)
arch_4=bestModel(Tubingen_arch[[4]])
##firstly,do PCA
pca_4=data.frame(arch_4$alphas)
colnames(pca_4)=c('A','B','C','D')
rownames(pca_4)=colnames(tübingen_cohort)
###working with degree of membership
pca_4=pca_4 %>%
  select(A,B,C,D) %>%
  mutate(
    arch= case_when(
      A > 0.6 ~ "archetypes A",
      B > 0.6 ~ "archetypes B",
      C > 0.6 ~ "archetypes C",
      D > 0.6 ~ "archetypes D",
      TRUE ~ 'mixed arch'
    )
  )

######do with ggplot
pca_tub=prcomp(pca_4[,-5])
fviz1<-fviz_pca_ind(pca_tub,
                     repel = TRUE,geom = "point")+geom_point(aes(colour=pca_input$arch,shape=pca_input$arch),size=3)+scale_color_brewer(palette="Dark2")+theme(
                       # Hide panel borders and remove grid lines
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()
                       )+
    theme(legend.text = element_text(size=15),legend.title = element_text(size=15),title = element_text(size=20))+
    labs(col= "Archetypes",shape="Archetypes")

pdf('pca_4arch.pdf')
fviz
dev.off()

nrow(pca_4[which(pca_4$arch=='archetypes A'),])
nrow(pca_4[which(pca_4$arch=='archetypes B'),])
nrow(pca_4[which(pca_4$arch=='archetypes C'),])
nrow(pca_4[which(pca_4$arch=='archetypes D'),])
nrow(pca_4[which(pca_4$arch=='mixed arch'),])

########pie chart
(dims = dim(pca_4))
results = pca_4%>%
  group_by(arch) %>%  
  summarize(Frequency = n()) %>% # return counts / frequencies
  mutate(Percent = paste0(round(Frequency / 72 * 100, 2), "%"))

pie=ggplot(results, aes(x = "", y = Frequency, fill=arch)) +
  geom_bar(width = 1, stat = "identity") + # this plots a stacked bar chart
  coord_polar(theta = "y", start = 0)+theme_void()+
  geom_text(aes(label = Percent),position = position_stack(vjust = 0.5))+scale_fill_brewer(palette="Dark2")

pdf('pie_4arch.pdf')
pie
dev.off()

#--------hratmap
#working with phenotypes 

pheno=fread('Thuebingen_Potsdam_cohorts_2outlier_removed.csv',skip = 1,data.table = FALSE)

pheno_tubingen=pheno %>% filter(str_detect(Sample_Plate,'A4726'))
pheno_short=pheno_tubingen %>% select(-Sample_Well,-Sample_Plate,-Pool_ID,-Sex,-Sample_Group,-Sentrix_Position,-Sentrix_ID,-cluster_tueftulip)

#### first, they must be  normally transformed 
pheno_short$age<-min_max_norm(pheno_short$age)
pheno_short$BMI<-min_max_norm(pheno_short$BMI)
pheno_short$WAIST<-min_max_norm(pheno_short$WAIST)
pheno_short$HIP<-min_max_norm(pheno_short$HIP)
pheno_short$CRP_mg_l<-min_max_norm(pheno_short$CRP_mg_l)
pheno_short$HDL_mmol_l<-min_max_norm(pheno_short$HDL_mmol_l)
pheno_short$Cholesterol_mmol_l<-min_max_norm(pheno_short$Cholesterol_mmol_l)
pheno_short$LDL_mmol_l<-min_max_norm(pheno_short$LDL_mmol_l)
pheno_short$GOT_AST_U_l=min_max_norm(pheno_short$GOT_AST_U_l)
pheno_short$GGT_U_l<-min_max_norm(pheno_short$GGT_U_l)
pheno_short$Creatinine_µmol_l<-min_max_norm(pheno_short$Creatinine_µmol_l)
pheno_short$MR_LFadj<-min_max_norm(pheno_short$MR_LFadj)
pheno_short$MR_SCAT_l<-min_max_norm(pheno_short$MR_SCAT_l)
pheno_short$MR_VAT_l<-min_max_norm(pheno_short$MR_VAT_l)
pheno_short$DI_calc<-min_max_norm(pheno_short$DI_calc)
pheno_short$ISI_calc<-min_max_norm(pheno_short$ISI_calc)



###filtering based one group 
#A=SAMPLE  Name
A=rownames(pca_4[which(pca_4$arch=='archetypes A'),])
B=rownames(pca_4[which(pca_4$arch=='archetypes B'),])
C=rownames(pca_4[which(pca_4$arch=='archetypes C'),])
D=rownames(pca_4[which(pca_4$arch=='archetypes D'),])
mixed=rownames(pca_4[which(pca_4$arch=='mixed arch'),])


####creating new data frame
#age
age_A=pheno_short$age[which(pheno_short$Sample_Name %in% A)]
age_B=pheno_short$age[which(pheno_short$Sample_Name %in% B)]
age_C=pheno_short$age[which(pheno_short$Sample_Name %in% C)]
age_D=pheno_short$age[which(pheno_short$Sample_Name %in% D)]
age_mixed=pheno_short$age[which(pheno_short$Sample_Name %in% mixed)]

library(coin)

max_length <- max(c(length(age_A), length(age_B),length(age_C),length(age_D),length(age_mixed)))    # Find out maximum length
max_length  
data <- data.frame(groupA= c(age_A,rep(NA, max_length - length(A))),groupB = c(age_B,rep(NA, max_length - length(B))),groupC= c(age_C,rep(NA, max_length - length(C))),groupD= c(age_D,rep(NA, max_length - length(D))),group_mixed= c(age_mixed,rep(NA, max_length - length(mixed))))


### pairwise comparison 
### group A :B,C,D,mixed
pvalue_A=wilcox.exact(data[,1],unlist(data[2:5]))$p.value
#group B : C,D,mixed
pvalue_B=wilcox.exact(data[,2],as.numeric(unlist(data[3:5])))$p.value
# group C : D,mixed
pvalue_C=wilcox.exact(data[,3],as.numeric(unlist(data[4:5])))$p.value
#Group D: mixed
pvalue_D=wilcox.exact(data[,4],as.numeric(unlist(data[,5])))$p.value

