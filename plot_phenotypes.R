library(ggbiplot)
library(ggpubr)

load('Load.Rdata')
data1=myLoad$pd
#pre
myLoad$pd$SEX='Female'
myLoad$pd$SEX=as.factor(myLoad$pd$SEX)
myLoad$pd$Saample_Name=as.factor(myLoad$pd$Sample_Name)
myLoad$pd$cluster_tueftulip=as.factor(myLoad$pd$cluster_tueftulip)
myLoad$pd$Slide=as.factor(myLoad$pd$Slide)
myLoad$pd=myLoad$pd[,-1]

#plot
p1=ggplot(data1,aes(x=cluster_tueftulip,y=WAIST,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
labs(x="Cluster",y= "WAIST")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p2=ggplot(data1,aes(x=cluster_tueftulip,y= HIP,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "HIP")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p3=ggplot(data1,aes(x=cluster_tueftulip,y= CRP_.mg.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "CRP_.mg.l")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p4=ggplot(data1,aes(x=cluster_tueftulip,y= HDL_.mmol.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "HDL_.mmol.l")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")


p5=ggplot(data1,aes(x=cluster_tueftulip,y=LDL_.mmol.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "LDL_.mmol.l.")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p6=ggplot(data1,aes(x=cluster_tueftulip,y= Cholesterol_.mmol.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "Cholesterol_.mmol.l.")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p7=ggplot(data1,aes(x=cluster_tueftulip,y= GOT.AST._.U.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "GOT.AST._.U.l.")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p8=ggplot(data1,aes(x=cluster_tueftulip,y= PRS,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "PRS")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

p9=ggplot(data1,aes(x=cluster_tueftulip,y= ISI.calc,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "ISI.calc")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9)




p11=ggplot(data1,aes(x=cluster_tueftulip,y= GGT_.U.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "GGT_.U.l.")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")


p22=ggplot(data1,aes(x=cluster_tueftulip,y= Creatinine_.µmol.l.,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "Creatinine_.µmol.l.")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")


p33=ggplot(data1,aes(x=cluster_tueftulip,y= MR_LFadj,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "MR_LFadj")+theme_classic(base_size = 25)+scale_color_discrete(name = 'Cluster')

p44=ggplot(data1,aes(x=cluster_tueftulip,y= MR_SCAT.l,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "MR_SCAT.l")+theme_classic(base_size = 25)+scale_color_discrete(name = 'Cluster')



p55=ggplot(data1,aes(x=cluster_tueftulip,y= MR_VAT.l,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "MR_VAT.l")+theme_classic(base_size = 25)+scale_color_discrete(name = 'Cluster')



p66=ggplot(data1,aes(x=cluster_tueftulip,y= DI.calc,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "DI.calc")+theme_classic(base_size = 25)+scale_color_discrete(name = 'Cluster')


p77=ggplot(data1,aes(x=cluster_tueftulip,y= age,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "Age")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")


p88=ggplot(data1,aes(x=cluster_tueftulip,y= BMI,color=cluster_tueftulip))+
  geom_boxplot()+geom_point(position = position_jitterdodge(), alpha=0.3)+
  labs(x="Cluster",y= "BMI")+theme_classic(base_size = 25)+scale_color_discrete(name = "Cluster")

ggarrange(p11,p22,p33,p44,p55,p66,p77,p88)












