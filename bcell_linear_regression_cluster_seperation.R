####linear regression


# libraries ----------------------------------------------------------------
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidylog)
library(janitor) #clean_names(df) or %>% clean_names()
library(ragg) #use agg_png to have quicker plots
library(report) #for session report
library(ggplot2)
require(svMisc)
library(tibble)
library(ggpubr)
# options -----------------------------------------------------------------
setwd('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction')
cell_com_tubingen=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction/cellcomposition.tsv')

phenodata <- fread("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/champ-3-outliers_females/raw/pheno.csv",skip=1)

phenodata=mutate(phenodata,Bcell=cell_com_tubingen$Bcell)
pheno_cleaned=phenodata   %>% 
  janitor::clean_names(
    case = "snake",
    replace = c("µ" = "micro", `'` = "", `"` = "", `%` = "_percent_", `#` = "_number_")) %>%
  mutate(cluster_whitehall = as.numeric(cluster_whitehall)) %>%
  mutate(cluster_tueftulip = as.numeric(cluster_tueftulip)) %>%
  mutate(waist = as.numeric(waist)) %>%
  mutate(hip = as.numeric(hip)) %>%
  mutate(mr_tat_l = as.numeric(mr_tat_l)) %>%
  mutate(glucose_000 = as.numeric(glucose_000)) %>%
  mutate(glucose_120 = as.numeric(glucose_120)) %>%
  mutate(glucose_auc = as.numeric(glucose_auc)) %>%
  mutate(cholesterol_mmol_l = as.numeric(cholesterol_mmol_l)) %>%
  mutate(ldl_mmol_l = as.numeric(ldl_mmol_l)) %>%
  mutate(triglycerides_mmol_l = as.numeric(triglycerides_mmol_l)) %>%
  mutate(got_ast_u_l = as.numeric(got_ast_u_l)) %>%
  mutate(gpt_alt_u_l = as.numeric(gpt_alt_u_l)) %>%
  mutate(ggt_u_l = as.numeric(ggt_u_l)) %>%
  mutate(creatinine_micromol_l = as.numeric(creatinine_micromol_l)) %>%
  mutate(igi_calc = as.numeric(igi_calc)) %>%
  dplyr::select(1:28,prs,bcell,-sex,-study_date,-sample_name,-cluster_whitehall,-age)


##take out na values
pheno_cleaned$mr_vat_l[is.na(pheno_cleaned$mr_vat_l)] <- mean(pheno_cleaned$mr_vat_l,na.rm = TRUE)
pheno_cleaned$mr_scat_l[is.na(pheno_cleaned$mr_scat_l)] <- mean(pheno_cleaned$mr_scat_l,na.rm = TRUE)
pheno_cleaned$mr_tat_l[is.na(pheno_cleaned$mr_tat_l)] <- mean(pheno_cleaned$mr_tat_l,na.rm = TRUE)




pheno_cleaned$cluster_tueftulip=as.character(pheno_cleaned$cluster_tueftulip)


###look at liver fat

###cluster sepration=

pheno_3=pheno_cleaned[pheno_cleaned$cluster_tueftulip=='3',]
pheno_6=pheno_cleaned[pheno_cleaned$cluster_tueftulip=='6',]
pheno_5=pheno_cleaned[pheno_cleaned$cluster_tueftulip=='5',]

fit5=lm(bcell ~ mr_l_fadj,data=pheno_5)

fit3=lm(bcell ~ mr_l_fadj,data=pheno_3 )

fit6=lm(bcell ~ mr_l_fadj,data=pheno_6 )

plt_5=ggplot(fit5$model, aes_string(x = names(fit5$model)[2], y = names(fit5$model)[1])) +
  labs(x = "liver fat") +
  labs(y = "Bcell") + 
  stat_smooth(method = "lm", col = "red", formula = "y ~ x") +
  scale_color_manual(values = c("#2aeb5b", "#ea5167", "#2dadcb", "#305a87")) +
  theme_light() +
  theme(plot.title = element_text(size = 10)) +labs(title = paste(
    "R2 = ", round(summary(fit5)$r.squared,3),
    " ,P =", round(summary(fit5)$coef[2, 4], 3),
    ',Cluster=',5
  ))+theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size=14),plot.title = element_text(size=14))


plt_3=ggplot(fit3$model, aes_string(x = names(fit3$model)[2], y = names(fit3$model)[1])) +
  labs(x = "liver fat") +
  labs(y = "Bcell") + 
  stat_smooth(method = "lm", col = "red", formula = "y ~ x") +
  scale_color_manual(values = c("#2aeb5b", "#ea5167", "#2dadcb", "#305a87")) +
  theme_light() +
  theme(plot.title = element_text(size = 10)) +labs(title = paste(
    "R2 = ", round(summary(fit3)$r.squared,3),
    " ,P =", round(summary(fit3)$coef[2, 4], 3),
    ',Cluster=',3
  ))+theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size=14),plot.title = element_text(size=14))


plt_6=ggplot(fit6$model, aes_string(x = names(fit6$model)[2], y = names(fit6$model)[1])) +
  labs(x = "liver fat") +
  labs(y = "Bcell") + 
  stat_smooth(method = "lm", col = "red", formula = "y ~ x") +
  scale_color_manual(values = c("#2aeb5b", "#ea5167", "#2dadcb", "#305a87")) +
  theme_light() +
  theme(plot.title = element_text(size = 10)) +labs(title = paste(
    "R2 = ", round(summary(fit6)$r.squared,3),
    " ,P =", round(summary(fit6)$coef[2, 4], 3),
    ',Cluster=',6
  ))+theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size=14),plot.title = element_text(size=14))


ggarrange(plt_3,plt_5,plt_6,ncol=3)
ggsave('liver_fat_clusters.pdf',width = 16,height =9,ggarrange(plt_3,plt_5,plt_6,ncol=3))






























