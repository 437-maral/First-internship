library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(factoextra)
library(readxl)
library(dplyr)
library(tibble)

load('myCombat.RData')

# first run dmp.R file 
#dc= discovery cohorts
# Rc= repliation cohorts

####dc
cpgs_intersect_35_36 <- intersect(my_combat_dmps_3_5$Name,my_combat_dmps_3_6$Name)#3697
cpgs_intersect_35_56 <- intersect(my_combat_dmps_3_5$Name,my_combat_dmps_6_5$Name)#6091
cpgs_intersect_36_56 <- intersect(my_combat_dmps_3_6$Name,my_combat_dmps_6_5$Name)



#####Rc

dmp3_5_dmp5_6=intersect(dmp3_5$Name,dmp5_6$Name)#19126
dmp5_6_dmp3_6=intersect(dmp5_6$Name,dmp3_6$Name)#28332
dmp3_6_dmp3_5=intersect(dmp3_5$Name,dmp3_6$Name)

### doing pairwise comparison 



pair_35_56=intersect(dmp3_5_dmp5_6,cpgs_intersect_35_56)##176
pair_35_36=intersect(cpgs_intersect_35_36,dmp3_6_dmp3_5)##161
pair_56_36=intersect(cpgs_intersect_36_56,dmp5_6_dmp3_6)#108


# 3,5 
pair_35=intersect(dmp3_5$Name,my_combat_dmps_3_5$Name)#4185
#5,6
pair_56=intersect(dmp5_6$Name,my_combat_dmps_6_5$Name)#2931
# 3,6
pair_36=intersect(dmp3_6$Name,my_combat_dmps_3_6$Name)# 2742

##new
com=as.data.frame(myCombat_array)
new_35=com[rownames(com) %in% pair_35,]
new_36=com[rownames(com) %in% pair_36,]
new_56=com[rownames(com) %in% pair_56,]

new35_56=com[rownames(com) %in% pair_35_56,]
new35_36=com[rownames(com) %in% pair_35_36,]
new56_36=com[rownames(com) %in% pair_56_36,]

### look for genes 
#3,5
anno_35=RatioSet(new_35, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes_35=data.frame(getAnnotation(anno_35))
left_35 <- left_join(rownames_to_column(new_35),rownames_to_column(probes_35))
genes_unique_3_5 <- left_35 %>% filter( UCSC_RefGene_Name != "") %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%  distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% pull(UCSC_RefGene_Name)%>% as.data.frame()

###imp
#### check the length for pathway
dim(genes_unique_3_5)
#[1] 2798    1

##3,6
anno_36=RatioSet(new_36, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes_36=data.frame(getAnnotation(anno_36))
left_36 <- left_join(rownames_to_column(new_36),rownames_to_column(probes_36))
genes_unique_3_6 <- left_36 %>% filter( UCSC_RefGene_Name != "") %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%  distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% pull(UCSC_RefGene_Name)%>% as.data.frame()

#####imp
#check
dim(genes_unique_3_6)#!!!!!!
#[1] 1895    1

#5,6
anno_56=RatioSet(new_56, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes_56=data.frame(getAnnotation(anno_56))
left_56 <- left_join(rownames_to_column(new_56),rownames_to_column(probes_56))
genes_unique_5_6 <- left_56%>% filter( UCSC_RefGene_Name != "") %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%  distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% pull(UCSC_RefGene_Name)%>% as.data.frame()

#######imp
# check 

dim(genes_unique_5_6)
#1967    1


anno56_36=RatioSet(new56_36, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes56_36=data.frame(getAnnotation(anno56_36))
left56_36 <- left_join(rownames_to_column(new56_36),rownames_to_column(probes56_36))
genes_unique56_36 <- left56_36 %>% filter( UCSC_RefGene_Name != "") %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%  distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% pull(UCSC_RefGene_Name)%>% as.data.frame()


dim(genes_unique56_36)
#77  1

anno56_35=RatioSet(new35_56, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probes56_35=data.frame(getAnnotation(anno56_35))
left56_35 <- left_join(rownames_to_column(new35_56),rownames_to_column(probes56_35))
genes_unique56_35 <- left56_35%>% filter( UCSC_RefGene_Name != "") %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%  distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% pull(UCSC_RefGene_Name)%>% as.data.frame()


dim(genes_unique56_35)
#140   1

anno36_35=RatioSet(new35_36, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
probe36_35=data.frame(getAnnotation(anno36_35))
left36_35 <- left_join(rownames_to_column(new35_36),rownames_to_column(probe36_35))
genes_unique36_35 <- left36_35%>% filter( UCSC_RefGene_Name != "") %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%  distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% pull(UCSC_RefGene_Name)%>% as.data.frame()





####
dim(genes_unique36_35)
#134   1

#### most important genes

write.table(x=genes_unique_3_5,file='genes_unique_3_5.tsv',sep='\t',row.names=F)
write.table(x=genes_unique_3_6,file='genes_unique_3_6.tsv',sep='\t',row.names=F)
write.table(x=genes_unique_5_6,file='genes_unique_5_6.tsv',sep='\t',row.names=F)


write.table(x=genes_unique36_35,file='genes_unique_36_35.tsv',sep='\t',row.names=F)
write.table(x=genes_unique56_35,file='genes_unique56_35.tsv',sep='\t',row.names=F)
write.table(x=genes_unique56_36,file='genes_unique56_36.tsv',sep='\t',row.names=F)






