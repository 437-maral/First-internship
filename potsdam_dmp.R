
library(ChAMP)
library(ggVennDiagram)
library(ggplot2)
setwd("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/Charlotte_Ling_method")

load('Load1.Rdata')#replication_cohort

load("new_combat.RData")



myLoad1$pd$new_3 <- gsub("5", "6", myLoad1$pd$cluster_tueftulip)
myLoad1$pd$new_5 <- gsub("3", "6", myLoad1$pd$cluster_tueftulip)
myLoad1$pd$new_6 <- gsub("5", "3", myLoad1$pd$cluster_tueftulip)

## Calculate DMPs
dmp_uni_3 <- champ.DMP(beta = new_mycombat,
                       pheno = as.character(myLoad1$pd$new_3),
                       compare.group = NULL,
                       adjPVal = 1,
                       adjust.method = "BH",
                       arraytype = "EPIC")

dmp_uni_5 <- champ.DMP(beta = new_mycombat,
                       pheno = as.character(myLoad1$pd$new_5),
                       compare.group = NULL,
                       adjPVal = 1,
                       adjust.method = "BH",
                       arraytype = "EPIC")

dmp_uni_6 <- champ.DMP(beta =new_mycombat,
                       pheno = as.character(myLoad1$pd$new_6),
                       compare.group = NULL,
                       adjPVal = 1,
                       adjust.method = "BH",
                       arraytype = "EPIC")



load('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/two_cohorts/Charlotte_Ling_method/potsdam_dmp.Rdata')

dmp3=dmp_uni_3$`6_to_3`[dmp_uni_3$`6_to_3`$P.Value< 0.05, ]
dmp5=dmp_uni_5$`6_to_5`[dmp_uni_5$`6_to_5`$P.Value< 0.05, ]
dmp6=dmp_uni_6$`6_to_3`[dmp_uni_6$`6_to_3`$P.Value< 0.05, ]


dim(dmp3[dmp3$adj.P.Val<0.05,])#5,503
dim(dmp5[dmp5$adj.P.Val<0.05,])#791
dim(dmp6[dmp6$adj.P.Val<0.05,])#2,198

dim(dmp3[dmp3$adj.P.Val<0.1,])#16,220
dim(dmp5[dmp5$adj.P.Val<0.1,])#3,777
dim(dmp6[dmp6$adj.P.Val<0.1,])#9,367

dim(dmp3[dmp3$adj.P.Val<0.2,])#47,220
dim(dmp5[dmp5$adj.P.Val<0.2,])#17,225
dim(dmp6[dmp6$adj.P.Val<0.2,])#34,681


dmp3$Name=rownames(dmp3)
dmp5$Name=rownames(dmp5)
dmp6$Name=rownames(dmp6)

length(unique(c(dmp3$Name,dmp6$Name,dmp5$Name)))

M=list('3'=dmp3$Name,'5'=dmp5$Name,'6'=dmp6$Name)

ggVennDiagram(M)


dim(dmp3)#5,503
dim(dmp5)#791
dim(dmp6)#2,198


###venn_diagram for different p_value

##0.05
dmp3_adj1=dmp3[dmp3$adj.P.Val<0.05,]#5,503

dmp5_adj1=dmp5[dmp5$adj.P.Val<0.05,]#791

dmp6_adj1=dmp6[dmp6$adj.P.Val<0.05,]#2,198

dmp3_adj1$Name=rownames(dmp3_adj1)
dmp5_adj1$Name=rownames(dmp5_adj1)
dmp6_adj1$Name=rownames(dmp6_adj1)



adj1=list('3'=dmp3_adj1$Name,'5'=dmp5_adj1$Name,'6'=dmp6_adj1$Name)


length(unique(c(dmp3_adj1$Name,dmp6_adj1$Name,dmp5_adj1$Name)))

ggVennDiagram(adj1)

##0.1

dmp3_adj2=dmp3[dmp3$adj.P.Val<0.1,]#5,503

dmp5_adj2=dmp5[dmp5$adj.P.Val<0.1,]#791

dmp6_adj2=dmp6[dmp6$adj.P.Val<0.1,]#2,198

dmp3_adj2$Name=rownames(dmp3_adj2)
dmp5_adj2$Name=rownames(dmp5_adj2)
dmp6_adj2$Name=rownames(dmp6_adj2)



adj2=list('3'=dmp3_adj2$Name,'5'=dmp5_adj2$Name,'6'=dmp6_adj2$Name)

ggVennDiagram(adj2)

length(unique(c(dmp3_adj2$Name,dmp6_adj2$Name,dmp5_adj2$Name)))
#0.2

dmp3_adj3=dmp3[dmp3$adj.P.Val<0.2,]#5,503

dmp5_adj3=dmp5[dmp5$adj.P.Val<0.2,]#791

dmp6_adj3=dmp6[dmp6$adj.P.Val<0.2,]#2,198

dmp3_adj3$Name=rownames(dmp3_adj3)
dmp5_adj3$Name=rownames(dmp5_adj3)
dmp6_adj3$Name=rownames(dmp6_adj3)



adj3=list('3'=dmp3_adj3$Name,'5'=dmp5_adj3$Name,'6'=dmp6_adj3$Name)

length(unique(c(dmp3_adj3$Name,dmp6_adj3$Name,dmp5_adj3$Name)))

ggVennDiagram(adj3)

pdf('venndiagaram_dmps.pdf')
ggVennDiagram(adj1)+ggtitle("p.value=0.05")+theme(plot.title=element_text(margin=margin(t=40,b=-30)))
ggVennDiagram(adj2)+ggtitle("p.value=0.1")+theme(plot.title=element_text(margin=margin(t=40,b=-30)))
ggVennDiagram(adj3)+ggtitle("p.value=0.2")+theme(plot.title=element_text(margin=margin(t=40,b=-30)))
dev.off()




