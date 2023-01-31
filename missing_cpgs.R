
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(factoextra)
library(readxl)
library(dplyr)
library(tibble)
library(openxlsx)


load('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/new_combat.Rdata')
load("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/!FINAL-results/2020CpgS.RData")
com=data.frame(new_mycombat)
`%!in%` <- Negate(`%in%`)
cpgs=my_combat_dmp_high_frq[rownames(my_combat_dmp_high_frq) %!in% rownames(com),]
cpgs_1980=data.frame(rownames(cpgs))


anno=RatioSet(cpgs, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes=data.frame(getAnnotation(anno))
left_35 <- left_join(rownames_to_column(cpgs),rownames_to_column(probes))
genes_unique_3_5 <- left_35%>% pull(UCSC_RefGene_Name)%>% as.data.frame()






load('two_combats.Rdata')



pheno1=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort//pheno.csv',skip = 1,data.table = FALSE)
pheno1$Sample_Name <- as.character(pheno1$Sample_Name)#like last time 
pheno1$cluster_tueftulip <- as.factor(pheno1$cluster_tueftulip)





Tübingen=two_myCombats[, colnames(two_myCombats) %in% pheno1$Sample_Name]

Tubingen2020=my_combat_dmp_high_frq[rownames(my_combat_dmp_high_frq) %!in% rownames(Tübingen),] #### 1998   72

cpgs_2020=data.frame(rownames(Tubingen2020))


anno_22=RatioSet(Tubingen2020, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes_22=data.frame(getAnnotation(anno_22))
left_22 <- left_join(rownames_to_column(Tubingen2020),rownames_to_column(probes_22))
genes_unique_22 <- left_22 %>% pull(UCSC_RefGene_Name)%>% as.data.frame()






