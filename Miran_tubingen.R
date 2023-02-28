library(data.table)
library(dplyr)
library(ggVennDiagram)
library(openxlsx)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(VennDiagram)
miRNA=fread('L:/!DIAB/Masoumeh Vojood/exp_validated_2023-02-22_T2D_subclusters.csv',data.table = TRUE)

#get mrian names
unique(miRNA$mir)
###we need ven diagram for these 4 genes miR-148a-3p
#[1] "hsa-let-7c-5p"   "hsa-miR-125b-5p" "hsa-miR-130a-3p" "hsa-miR-148a-3p" "hsa-miR-192-5p" 
#[6] "hsa-miR-19b-3p"  "hsa-miR-342-3p"  "hsa-miR-532-5p"  "hsa-miR-874-3p" 

###between clusters they are more significant
#miR-148a-3p
#miR-874-3p
#mir-192-5p
#mir-125b-5b

#####

##do filtering

#"hsa-let-7c-5p" 

miRNA_hsa_let=miRNA[miRNA$mir =='hsa-let-7c-5p',]#2300  
###without any duplication
select(miRNA_hsa_let,gene_ID) %>% unique() %>% dim()# 2298
select(miRNA_hsa_let,Symbol) %>% unique() %>% dim()# 2238 
###pull genes names
select(miRNA_hsa_let,Symbol) %>% unique() %>% pull() %>% data.frame()


#hsa-miR-125b-5p

miRNA_hsa_125=miRNA[miRNA$mir =='hsa-miR-125b-5p',]# 1421
###without any duplication
select(miRNA_hsa_125,gene_ID) %>% unique() %>% dim()# 1420 
select(miRNA_hsa_125,Symbol) %>% unique() %>% dim()# 1376
###pull genes names
select(miRNA_hsa_125,Symbol) %>% unique() %>% pull() %>% data.frame()

#'hsa-miR-130a-3p'

miRNA_hsa_130=miRNA[miRNA$mir =='hsa-miR-130a-3p',]#  2126  
###without any duplication
select(miRNA_hsa_130,gene_ID) %>% unique() %>% dim()# 2124 
select(miRNA_hsa_130,Symbol) %>% unique() %>% dim()# 2099 
###pull genes names
select(miRNA_hsa_130,Symbol) %>% unique() %>% pull() %>% data.frame()

#hsa-miR-148a-3p"

miRNA_hsa_148=miRNA[miRNA$mir =='hsa-miR-148a-3p',]#  1576 
###without any duplication
select(miRNA_hsa_148,gene_ID) %>% unique() %>% dim()# 1576 
select(miRNA_hsa_148,Symbol) %>% unique() %>% dim()# 1545 
###pull genes names
select(miRNA_hsa_148,Symbol) %>% unique() %>% pull() %>% data.frame()



#hsa-miR-192-5p"

miRNA_hsa_192=miRNA[miRNA$mir =='hsa-miR-192-5p',]# 1625 
###without any duplication
select(miRNA_hsa_192,gene_ID) %>% unique() %>% dim()# 1597
select(miRNA_hsa_192,Symbol) %>% unique() %>% dim()# 1396
###pull genes names
select(miRNA_hsa_192,Symbol) %>% unique() %>% pull() %>% data.frame()



#hsa-miR-19b-3p"

miRNA_hsa_19=miRNA[miRNA$mir =='hsa-miR-19b-3p',]# 1991 
###without any duplication
select(miRNA_hsa_19,gene_ID) %>% unique() %>% dim()#  1990
select(miRNA_hsa_19,Symbol) %>% unique() %>% dim()# 1955
###pull genes names

select(miRNA_hsa_532,Symbol) %>% unique() %>% pull() %>% data.frame()

#hsa-miR-342-3p"

miRNA_hsa_342=miRNA[miRNA$mir =='hsa-miR-342-3p',]# 866 
###without any duplication
select(miRNA_hsa_342,gene_ID) %>% unique() %>% dim()#  866   
select(miRNA_hsa_342,Symbol) %>% unique() %>% dim()# 825 
###pull genes names
select(miRNA_hsa_342,Symbol) %>% unique() %>% pull()%>% data.frame()

#hsa-miR-532-5p"

miRNA_hsa_532=miRNA[miRNA$mir =='hsa-miR-532-5p',]# 372  
###without any duplication
select(miRNA_hsa_532,gene_ID) %>% unique() %>% dim()#  372   
select(miRNA_hsa_532,Symbol) %>% unique() %>% dim()#371 
###pull genes names
select(miRNA_hsa_532,Symbol) %>% unique() %>% pull() %>% data.frame()


"hsa-miR-874-3p" 

miRNA_hsa_874=miRNA[miRNA$mir =="hsa-miR-874-3p" ,]# 370
###without any duplication
select(miRNA_hsa_874,gene_ID) %>% unique() %>% dim()#  370
select(miRNA_hsa_874,Symbol) %>% unique() %>% dim()#368
select(miRNA_hsa_874,Symbol) %>% unique() %>% pull() %>% data.frame()


##do overlapiing for four miRNA

#miR-148a-3p n=1576 
#miR-874-3p n=370
#mir-192-5p n=1625 
#mir-125b-5b n=1421


###do overlapping for four mirna




m=list('miR-148a-3p'=miRNA_hsa_148$Symbol,'miR-874-3p'=miRNA_hsa_874$Symbol,'mir-192-5p'=miRNA_hsa_192$Symbol,'mir-125b-5b'=miRNA_hsa_125$Symbol)

plt=ggVennDiagram(m)+ggplot2::scale_fill_gradient(low="white",high = "white")

###do overlapping for four mirna

overlapping_mirna<- Reduce(intersect,list(miRNA_hsa_148$Symbol,miRNA_hsa_874$Symbol,miRNA_hsa_192$Symbol,miRNA_hsa_125$Symbol))




