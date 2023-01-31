
library(ggVennDiagram)
library(data.table)


##directory

setwd('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/two_cohorts')
#function
source("L:/!DIAB/Masoumeh Vojood/functions/PAM.R")
source("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/codes/dife_pca_2021_changes.R")


###loading
pheno2<-fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno.csv",skip=1,data.table=FALSE)#here the 1 outlier is still in 45086 
pheno1=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort//pheno.csv',skip = 1,data.table = FALSE)
load('two_combats_2outliersrm.Rdata')


##potsdam

load('dmp_twocohort_potsdam_2outliersrm.Rdata')

dmp6_5=dmptest_potsdam$`6_to_5`[dmptest_potsdam$`6_to_5`$P.Value<0.05,]
dmp3_6=dmptest_potsdam$`6_to_3`[dmptest_potsdam$`6_to_3`$P.Value<0.05,]
dmp5_3=dmptest_potsdam$`5_to_3`[dmptest_potsdam$`5_to_3`$P.Value<0.05,]

##add Name column for cpgs
dmp6_5$Name=rownames(dmp6_5)
dmp3_6$Name=rownames(dmp3_6)
dmp5_3$Name=rownames(dmp5_3)


all_cpgs_po=unique(c(dmp5_3$Name,dmp3_6$Name,dmp6_5$Name))#85603



M <-list('3 vs. 5'=dmp5_3$Name,
         '3 vs. 6'=dmp3_6$Name,
         '5 vs. 6'=dmp6_5$Name)
ggVennDiagram(M,color = 1, lwd = 0.8, lty = 1)



###pam 
##all dmps



Potsdam_cohort=two_myCombats[,colnames(two_myCombats) %in% pheno2$Sample_Name]


all_cpgs_po=unique(c(dmp5_3$Name,dmp3_6$Name,dmp6_5$Name))#85603

mycombat_alldmps=Potsdam_cohort[rownames(Potsdam_cohort) %in% all_cpgs_po,]###85603    46
mycombat_alldmps=data.frame(mycombat_alldmps)


names3 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 6)])


colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names3)]<-paste0(colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names3)],"_3")
colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names5)]<-paste0(colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names5)],"_5")
colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names6)]<-paste0(colnames(mycombat_alldmps)[which(colnames(mycombat_alldmps) %in% names6)],"_6")
colnames(mycombat_alldmps)



pam_function(input=mycombat_alldmps,class="unknown4")#69.57"

####intersecting
dmps <- c(dmp6_5$Name,dmp5_3$Name,dmp3_6$Name)
duplicated  <- dmps[duplicated(dmps)] #22k



mycombat_intersecting=Potsdam_cohort[rownames(Potsdam_cohort) %in% duplicated ,]#20453    46

mycombat_intersecting=data.frame(mycombat_intersecting)


names3 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 6)])

colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names3)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names3)],"_3")
colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names5)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names5)],"_5")
colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names6)]<-paste0(colnames(mycombat_intersecting)[which(colnames(mycombat_intersecting) %in% names6)],"_6")
colnames(mycombat_intersecting)



pam_function(input=mycombat_intersecting,class='unknown2')#73.91"

###specific dmps
duplicated  <- dmps[duplicated(dmps)] #22
specific_dmps <- subset(dmps, !dmps %in% duplicated) #65k


mycombat_spec=Potsdam_cohort[rownames(Potsdam_cohort) %in% specific_dmps,]#  65150 
mycombat_spec=data.frame(mycombat_spec)




names3 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 6)])


colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names3)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names3)],"_3")
colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names5)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names5)],"_5")
colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names6)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names6)],"_6")
colnames(mycombat_spec)



pam_function(input=mycombat_spec,class="unknown2")###7


####overlapped
new=Reduce(intersect,list(dmp5_3$Name,dmp3_6$Name,dmp6_5$Name))


mycombat_overlapped=Potsdam_cohort[rownames(Potsdam_cohort) %in% new,]#  59
mycombat_overlapped=data.frame(mycombat_overlapped)


names3 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno2$Sample_Name[which(pheno2$cluster_tueftulip == 6)])


colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names3)]<-paste0(colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names3)],"_3")
colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names5)]<-paste0(colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names5)],"_5")
colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names6)]<-paste0(colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names6)],"_6")
colnames(mycombat_overlapped)




pam_function(input=mycombat_overlapped,class="unknown2")###7


#####pca
cluster_rep_alldmps=c(colnames(mycombat_alldmps))
cluster3_rep<- c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '3' )])
cluster5_rep<- c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '5')])
cluster6_rep=c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '6')])

######
cluster_rep_alldmps[which(colnames(mycombat_alldmps) %in% cluster3_rep)]<-'3'
cluster_rep_alldmps[which(colnames(mycombat_alldmps) %in% cluster5_rep)]<-'5'
cluster_rep_alldmps[which(colnames(mycombat_alldmps) %in% cluster6_rep)]<-'6'

PCA_rep_cluster=dife_pca(mycombat_alldmps,colorgroups = cluster_rep_alldmps)



l=biplot(PCA_rep_cluster,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Potsdam Cohort')

#####pca for intersecting dmps

cluster_rep_intersecting=c(colnames(mycombat_intersecting))
cluster3_rep2<- c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '3' )])
cluster5_rep2<- c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '5')])
cluster6_rep2=c(pheno2$Sample_Name[which(pheno2$cluster_tueftulip ==  '6')])

######
cluster_rep_intersecting[which(colnames(mycombat_intersecting) %in% cluster3_rep2)]<-'3'
cluster_rep_intersecting[which(colnames(mycombat_intersecting) %in% cluster5_rep2)]<-'5'
cluster_rep_intersecting[which(colnames(mycombat_intersecting) %in% cluster6_rep2)]<-'6'

PCA_rep_cluste_intersecting=dife_pca(mycombat_intersecting,colorgroups = cluster_rep_intersecting)

l2=biplot(PCA_rep_cluste_intersecting,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Potsdam Cohort')




#------------------------------------------------------------------
#####Tubinegn

load('dmp_twocohort_tubingen_2outliersrm.Rdata')

dmp6_5_tue=dmptest_tubingen$`6_to_5`[dmptest_tubingen$`6_to_5`$P.Value<0.05,]
dmp3_6_tue=dmptest_tubingen$`3_to_6`[dmptest_tubingen$`3_to_6`$P.Value<0.05,]
dmp3_5_tue=dmptest_tubingen$`3_to_5`[dmptest_tubingen$`3_to_5`$P.Value<0.05,]


##add Name column for cpgs
dmp6_5_tue$Name=rownames(dmp6_5_tue)
dmp3_6_tue$Name=rownames(dmp3_6_tue)
dmp3_5_tue$Name=rownames(dmp3_5_tue)

#####pam

Tubingen_cohort=two_myCombats[,colnames(two_myCombats) %in% pheno1$Sample_Name]

all_cpgs_tu=unique(c(dmp6_5_tue$Name,dmp3_6_tue$Name,dmp3_5_tue$Name))#65895

mycombat_alldmps_tue=Tubingen_cohort[rownames(Tubingen_cohort) %in% all_cpgs_tu,]#65895    72
mycombat_alldmps_tue=data.frame(mycombat_alldmps_tue)

names3_tue<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 3)])
names5_tue<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 5)])
names6_tue<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 6)])


colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names3_tue)]<-paste0(colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names3_tue)],"_3")
colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names5_tue)]<-paste0(colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names5_tue)],"_5")
colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names6_tue)]<-paste0(colnames(mycombat_alldmps_tue)[which(colnames(mycombat_alldmps_tue) %in% names6_tue)],"_6")
colnames(mycombat_alldmps_tue)




pam_function(input=mycombat_alldmps_tue,class="sofya")#52%


dmps_tue<- c(dmp6_5_tue$Name,dmp3_6_tue$Name,dmp3_5_tue$Name)
duplicated_tue<- dmps_tue[duplicated(dmps_tue)]#12901



mycombat_intersecting_tue=Tubingen_cohort[rownames(Tubingen_cohort) %in% duplicated_tue ,]#12885    72


mycombat_intersecting_tue=data.frame(mycombat_intersecting_tue)

names3_tue<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 3)])
names5_tue<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 5)])
names6_tue<- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 6)])


colnames(mycombat_intersecting_tue)[which(colnames(mycombat_intersecting_tue) %in% names3_tue)]<-paste0(colnames(mycombat_intersecting_tue)[which(colnames(mycombat_intersecting_tue) %in% names3_tue)],"_3")
colnames(mycombat_intersecting_tue)[which(colnames(mycombat_intersecting_tue) %in% names5_tue)]<-paste0(colnames(mycombat_intersecting_tue)[which(colnames(mycombat_intersecting_tue) %in% names5_tue)],"_5")
colnames(mycombat_intersecting_tue)[which(colnames(mycombat_intersecting_tue) %in% names6_tue)]<-paste0(colnames(mycombat_intersecting_tue)[which(colnames(mycombat_intersecting_tue) %in% names6_tue)],"_6")
colnames(mycombat_intersecting_tue)


pam_function(input=mycombat_intersecting_tue,class="markus")#66%


###specified
dmps_tue<- c(dmp6_5_tue$Name,dmp3_6_tue$Name,dmp3_5_tue$Name)
duplicated  <- dmps_tue[duplicated(dmps_tue)] #22
specific_dmps <- subset(dmps_tue, !dmps_tue %in% duplicated) #65k


mycombat_spec=Tubingen_cohort[rownames(Tubingen_cohort) %in% specific_dmps,]#  65150 
mycombat_spec=data.frame(mycombat_spec)




names3 <- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 6)])


colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names3)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names3)],"_3")
colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names5)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names5)],"_5")
colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names6)]<-paste0(colnames(mycombat_spec)[which(colnames(mycombat_spec) %in% names6)],"_6")
colnames(mycombat_spec)



pam_function(input=mycombat_spec,class="unknown4")###7



####overlapped region 


new=Reduce(intersect,list(dmp6_5_tue$Name,dmp3_6_tue$Name,dmp3_5_tue$Name))


mycombat_overlapped=Tubingen_cohort[rownames(Tubingen_cohort) %in% new,]#  65150 
mycombat_overlapped=data.frame(mycombat_overlapped)


names3 <- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 3)])
names5 <- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 5)])
names6 <- paste0("X",pheno1$Sample_Name[which(pheno1$cluster_tueftulip == 6)])


colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names3)]<-paste0(colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names3)],"_3")
colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names5)]<-paste0(colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names5)],"_5")
colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names6)]<-paste0(colnames(mycombat_overlapped)[which(colnames(mycombat_overlapped) %in% names6)],"_6")
colnames(mycombat_overlapped)




pam_function(input=mycombat_overlapped,class="unknown4")###7

###PCA

cluster_disc_tue=c(colnames(mycombat_alldmps_tue))
cluster3_disc<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '3' )])
cluster5_disc<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '5')])
cluster6_disc<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '6')])

######
cluster_disc_tue[which(colnames(mycombat_alldmps_tue) %in% cluster3_disc)]<-'3'
cluster_disc_tue[which(colnames(mycombat_alldmps_tue) %in% cluster5_disc)]<-'5'
cluster_disc_tue[which(colnames(mycombat_alldmps_tue) %in% cluster6_disc)]<-'6'

PCA_disc_cluster=dife_pca(mycombat_alldmps_tue,colorgroups = cluster_disc_tue)



f=biplot(PCA_disc_cluster,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Tubinegn')


cluster_disc_intersecting=c(colnames(mycombat_intersecting_tue))
cluster3_disc2<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '3' )])
cluster5_dis2c<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '5')])
cluster6_disc2<- c(pheno1$Sample_Name[which(pheno1$cluster_tueftulip ==  '6')])

######
cluster_disc_intersecting[which(colnames(mycombat_intersecting_tue) %in% cluster3_disc)]<-'3'
cluster_disc_intersecting[which(colnames(mycombat_intersecting_tue) %in% cluster5_disc)]<-'5'
cluster_disc_intersecting[which(colnames(mycombat_intersecting_tue) %in% cluster6_disc)]<-'6'

PCA_disc_intersecting=dife_pca(mycombat_intersecting_tue,colorgroups = cluster_disc_intersecting)

f2=biplot(PCA_disc_intersecting,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Tubinegn')


pdf('all_DMPs.pdf')
print(biplot(PCA_rep_cluster,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Potsdam Cohort'))
print(biplot(PCA_disc_cluster,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Tubinegn'))
dev.off()


pdf('intersecting_DMPs.pdf')
print(biplot(PCA_rep_cluste_intersecting,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Potsdam Cohort'))
print(biplot(PCA_disc_intersecting,colby ='group',hline = 0, vline = 0,legendPosition = 'right',title='Tubinegn'))
dev.off()






