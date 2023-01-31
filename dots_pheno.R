library(ggbiplot)
library(ggpubr)


pheno1=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort//pheno.csv',skip = 1,data.table = FALSE)
pheno1$Sample_Name <- as.character(pheno1$Sample_Name)#like last time 
pheno1$cluster_tueftulip <- as.factor(pheno1$cluster_tueftulip)

pheno2<-fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/replication_cohort/raw/pheno1.csv",skip = 1,data.table = FALSE)
pheno2$Sample_Name <- as.character(pheno2$Sample_Name)#like last time 
pheno2$cluster_tueftulip <- as.factor(pheno2$cluster_tueftulip)


pheno1$cohorts='Tübingen'
pheno2$cohorts='Potsdam'
pheno1=select(pheno1,- "cluster_whitehall",-"study_date" ,-"MR_TAT_l",- "glucose_000",-"glucose_120",-"glucose_AUC",-"Triglycerides_mmol_l",
              -"class_02arb" ,-"class_01acei" ,-"class_03thiazids" ,-"class_04otherdiuretics", -"Basename",-"GPT_ALT_U_l"  )

pheno1=select(pheno1,-"class_05betablockers",-"class_06statins")           
pheno1=select(pheno1,-"IGI_calc")
pheno1=dplyr::rename(pheno1,"CRP_[mg/l]" ="CRP_mg_l" ,"HDL_[mmol/l]" ="HDL_mmol_l", "LDL_[mmol/l]" = "LDL_mmol_l" ,      
                     "Cholesterol_[mmol/l]" ="Cholesterol_mmol_l","GOT(AST)_[U/l]" =  "GOT_AST_U_l" ,"GGT_[U/l]"="GGT_U_l" , 
                     "Creatinine_[µmol/l]" ="Creatinine_µmol_l",
                     "MR_LFadj"=  "MR_LFadj", "MR_SCAT/l"="MR_SCAT_l", "MR_VAT/l" =   "MR_VAT_l","DI.calc"= "DI_calc",         
                     "ISI.calc"="ISI_calc")
new=rbind(pheno2,pheno1)     



p1=ggplot(new,aes(x = cluster_tueftulip, y =ISI.calc,fill=factor(cohorts)))+geom_dotplot(binaxis ='y', dotsize = 1.15,stackdir = "center")+facet_wrap(vars(cohorts))

       
       
p2=ggplot(new,aes(x = cluster_tueftulip, y =PRS,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',dotsize = 1.15,stackdir = "center")+facet_wrap(vars(cohorts))



p3=ggplot(new,aes(x = cluster_tueftulip, y =BMI,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+theme_bw()+theme(axis.text = element_text(size = 18),axis.title=element_text(size = 18),plot.title = element_text(size = 18),legend.text = element_text(size=18),strip.background = element_rect(colour="black", fill="white"),legend.title = element_text(size=18),strip.text = element_text(size=12))+facet_wrap(vars(cohorts))+labs(fill = "Cohorts")


p4=ggplot(new,aes(x = cluster_tueftulip, y =`Creatinine_[µmol/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))



p5=ggplot(new,aes(x = cluster_tueftulip, y =`LDL_[mmol/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))


p6=ggplot(new,aes(x = cluster_tueftulip, y =`MR_LFadj`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))

p7=ggplot(new,aes(x = cluster_tueftulip, y =age,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))


p8=ggplot(new,aes(x = cluster_tueftulip, y =WAIST,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))

p9=ggplot(new,aes(x = cluster_tueftulip, y =HIP,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))

p10=ggplot(new,aes(x = cluster_tueftulip, y =`MR_SCAT/l`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))

p11=ggplot(new,aes(x = cluster_tueftulip, y =`MR_VAT/l`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))


p12=ggplot(new,aes(x = cluster_tueftulip, y =`CRP_[mg/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))

p13=ggplot(new,aes(x = cluster_tueftulip, y =`Cholesterol_[mmol/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))



p14=ggplot(new,aes(x = cluster_tueftulip, y =`HDL_[mmol/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))


p15=ggplot(new,aes(x = cluster_tueftulip, y =`GOT(AST)_[U/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))


p16=ggplot(new,aes(x = cluster_tueftulip, y =`GGT_[U/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))


p17=ggplot(new,aes(x = cluster_tueftulip, y =DI.calc,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))




ggarrange(p1,p2,p3,p4)

ggsave('plot1.png',width = 15, height = 8 ,dpi = 100)

ggarrange(p5,p6,p7,p8)

ggsave('plot2.png',width = 15, height = 8 ,dpi = 100)

ggarrange(p9,p10,p11,p12)

ggsave('plot3.png',width = 15, height = 8 ,dpi = 100)

ggarrange(p13,p14,p15)

ggsave('plot4.png',width = 15, height = 8 ,dpi = 100)

ggarrange(p16,p17,nrow=2,ncol=2)

ggsave('plot4.png',width = 15, height = 8 ,dpi = 100)


##### tp resent

p33=ggplot(new,aes(x = cluster_tueftulip, y =BMI,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+theme_bw()+theme(axis.text = element_text(size = 18),axis.title=element_text(size = 18),plot.title = element_text(size = 18),strip.background = element_rect(colour="black", fill="white"),legend.position="none",strip.text = element_text(size=12))+facet_wrap(vars(cohorts))

p44=ggplot(new,aes(x = cluster_tueftulip, y =`Creatinine_[µmol/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+theme_bw()+theme(axis.text = element_text(size = 18),axis.title=element_text(size = 18),plot.title = element_text(size = 18),strip.background = element_rect(colour="black"), fill="white",legend.position="none",strip.text = element_text(size=12))+facet_wrap(vars(cohorts))

p55=ggplot(new,aes(x = cluster_tueftulip, y =`GGT_[U/l]`,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+facet_wrap(vars(cohorts))+theme_bw()+theme(axis.text = element_text(size = 18),axis.title=element_text(size = 18),plot.title = element_text(size = 18),legend.text = element_text(size=18),strip.background = element_rect(colour="black", fill="white"),legend.title = element_text(size=18),strip.text = element_text(size=12))+facet_wrap(vars(cohorts))+labs(fill = "Cohorts")


ggarrange(p33,p44,p55,nrow=3,ncol=2)


p11=ggplot(new,aes(x = cluster_tueftulip, y =HIP,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+theme_bw()+theme(axis.text = element_text(size = 18),axis.title=element_text(size = 18),plot.title = element_text(size = 18),legend.text = element_text(size=18),strip.background = element_rect(colour="black", fill="white"),legend.title = element_text(size=18),strip.text = element_text(size=12))+facet_wrap(vars(cohorts))+labs(fill = "Cohorts")


p77=ggplot(new,aes(x = cluster_tueftulip, y =age,fill=factor(cohorts)))+geom_dotplot(binaxis ='y',stackdir = "center")+theme_bw()+theme(axis.text = element_text(size = 18),axis.title=element_text(size = 18),plot.title = element_text(size = 18),legend.text = element_text(size=18),strip.background = element_rect(colour="black", fill="white"),legend.title = element_text(size=18),strip.text = element_text(size=12))+facet_wrap(vars(cohorts))+labs(fill = "Cohorts")


ggarrange(p11,p77,nrow=3,ncol=2)
