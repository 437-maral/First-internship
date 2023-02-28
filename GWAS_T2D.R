##in this script I linked the SNPs with cpgs

###I liked the chromosome position form GWAS  data and SNP position FROM QTL data

##I found methyl QTL in this website 'https://epicmeqtl.kcl.ac.uk/'

setwd('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/gwas')
library(data.table)




##load data 
####start with type 2 diabates
diabetes_snp=read.xlsx('QTL_data.xlsx')
diabetes_gwas=read.xlsx('type2_diabetes.xlsx')

unique(diabetes_gwas$`DISEASE/TRAIT`)###different type of traits associated T2D , but I look at type 2 diabetes

t2d=diabetes_gwas[diabetes_gwas$`DISEASE/TRAIT` %in% c("Type 2 diabetes","Severe insulin-resistant type 2 diabetes","Severe insulin-deficient type 2 diabetes","Triglyceride levels in type 2 diabetes"),]

##match
cpgs_t2d=diabetes_snp$CpG[t2d$CHR_POS %in% diabetes_snp$SNP.pos]


###look at BMI 

BMI_gwas=read.xlsx('BMI.xlsx')

##look at BMI
unique(BMI_gwas$`DISEASE/TRAIT`)####I look at the interaction between T2D and VMI and BMI itself

BMI=BMI_gwas[BMI_gwas$`DISEASE/TRAIT` %in% c("Body mass index"  , "Body mass index and type 2 diabetes (pairwise)"),]

cpgs_BMI=diabetes_snp$CpG[BMI$CHR_POS %in% diabetes_snp$SNP.pos]


###venndiagram

library(ggVennDiagram)

ggVennDiagram(list('T2D'=cpgs_t2d,'BMI'=cpgs_BMI))









#####look at liver fat
#look at another datatset

liver_gwas=read.xlsx('liver_fat.xlsx')
unique(liver_gwas$`DISEASE/TRAIT`)###liver which disease 

liver=liver_gwas[liver_gwas$`DISEASE/TRAIT` %in% c("Percent liver fat","Liver fat" ,"Liver fat content (MRI proton density fat fraction measure)") ,]

cpgs_liver=diabetes_snp$CpG[liver_gwas$CHR_POS %in% diabetes_snp$SNP.pos]


##insulin_secreation


insulin_secretion_gwas=read.xlsx('insulin_secretion.xlsx')
unique(insulin_secretion_gwas$`DISEASE/TRAIT`)###insulin which disease 

inslulin_secretion=insulin_secretion_gwas[insulin_secretion_gwas$`DISEASE/TRAIT`=='Insulin secretion rate' ,]

cpgs_inslulin_secretion=diabetes_snp$CpG[insulin_secretion_gwas$CHR_POS %in% diabetes_snp$SNP.pos]


##insulin resistance 

insulin_resistance_gwas=read.xlsx('insulin_resistance.xlsx')
unique(insulin_resistance_gwas$`DISEASE/TRAIT`)###insulin which disease 

inslulin_resistance=insulin_resistance_gwas[insulin_resistance_gwas$`DISEASE/TRAIT`%in% c("Insulin resistance/response" ,"Diabetes related insulin traits" ,"Insulin-related traits" ) ,]

cpgs_inslulin_resistance=diabetes_snp$CpG[insulin_resistance_gwas$CHR_POS %in% diabetes_snp$SNP.pos]




