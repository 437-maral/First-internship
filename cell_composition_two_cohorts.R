

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

cell_com_tubingen=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction/cellcomposition.tsv')


cell_com_potsdam=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/cellcomposition.tsv')


my_comparisons <- list(c('Potsdam','Tübingen'))

cell_com_tubingen=select(cell_com_tubingen,-Name)
cell_com_tubingen$Group=as.character(cell_com_tubingen$Group)
cell_com_tubingen$cohort='Tübingen'

cell_com_potsdam=select(cell_com_potsdam,-Name)
cell_com_potsdam$Group=as.character(cell_com_potsdam$Group)
cell_com_potsdam$cohort='Potsdam'

combined_data=rbind(cell_com_tubingen,cell_com_potsdam)


my_comparisons <- list(c('3.Potsdam','3.Tübingen'),c('5.Potsdam','5.Tübingen'),c('6.Potsdam','6.Tübingen'),
                       c('3.Potsdam','5.Potsdam'),c('3.Potsdam','6.Potsdam'),c('3.Tübingen','5.Tübingen'),c('3.Tübingen','6.Tübingen'),
                       c('5.Potsdam','6.Potsdam'),c('5.Tübingen','6.Tübingen'))

cd8t=combined_data%>%
  ggplot(aes(x=interaction(Group,cohort),y=CD8T,fill=cohort)) + 
  geom_boxplot(aes(group = interaction(Group,cohort)), outlier.shape = NA,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
  scale_x_discrete('',limits=c('3.Potsdam','3.Tübingen','5.Potsdam','5.Tübingen','6.Potsdam','6.Tübingen'),
                   labels=c('3.P','3.T','5.P','5.T','6.P','6.T'))+
  stat_compare_means(comparisons=my_comparisons)+
  scale_fill_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle('CD8T')


CD4T=combined_data%>%
  ggplot(aes(x=interaction(Group,cohort),y=CD4T,fill=cohort)) + 
  geom_boxplot(aes(group = interaction(Group,cohort)), outlier.shape = NA,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
  scale_x_discrete('',limits=c('3.Potsdam','3.Tübingen','5.Potsdam','5.Tübingen','6.Potsdam','6.Tübingen'),
                   labels=c('3.P','3.T','5.P','5.T','6.P','6.T'))+
  stat_compare_means(comparisons=my_comparisons)+
  scale_fill_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle('CD4T')




NK=combined_data%>%
  ggplot(aes(x=interaction(Group,cohort),y=NK,fill=cohort)) + 
  geom_boxplot(aes(group = interaction(Group,cohort)), outlier.shape = NA,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
  scale_x_discrete('',limits=c('3.Potsdam','3.Tübingen','5.Potsdam','5.Tübingen','6.Potsdam','6.Tübingen'),
                   labels=c('3.P','3.T','5.P','5.T','6.P','6.T'))+
  stat_compare_means(comparisons=my_comparisons)+
  scale_fill_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle('NK')



Bcell=combined_data%>%
  ggplot(aes(x=interaction(Group,cohort),y=Bcell,fill=cohort)) + 
  geom_boxplot(aes(group = interaction(Group,cohort)), outlier.shape = NA,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
  scale_x_discrete('',limits=c('3.Potsdam','3.Tübingen','5.Potsdam','5.Tübingen','6.Potsdam','6.Tübingen'),
                   labels=c('3.P','3.T','5.P','5.T','6.P','6.T'))+
  stat_compare_means(comparisons=my_comparisons)+
  scale_fill_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle('Bcell')


Mono=combined_data%>%
  ggplot(aes(x=interaction(Group,cohort),y=Mono,fill=cohort)) + 
  geom_boxplot(aes(group = interaction(Group,cohort)), outlier.shape = NA,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
  scale_x_discrete('',limits=c('3.Potsdam','3.Tübingen','5.Potsdam','5.Tübingen','6.Potsdam','6.Tübingen'),
                   labels=c('3.P','3.T','5.P','5.T','6.P','6.T'))+
  stat_compare_means(comparisons=my_comparisons)+
  scale_fill_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle('Mono')



Gran=combined_data%>%
  ggplot(aes(x=interaction(Group,cohort),y=Gran,fill=cohort)) + 
  geom_boxplot(aes(group = interaction(Group,cohort)), outlier.shape = NA,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
  scale_x_discrete('',limits=c('3.Potsdam','3.Tübingen','5.Potsdam','5.Tübingen','6.Potsdam','6.Tübingen'),
                   labels=c('3.P','3.T','5.P','5.T','6.P','6.T'))+
  stat_compare_means(comparisons=my_comparisons)+
  scale_fill_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle('Gran')




l=ggarrange(cd8t,CD4T,Bcell,NK,Mono,Gran,nrow = 2, ncol = 3)
ggsave(filename = 'L:/Bioinformatics_Unit/T2D_Subcluster_Studie/tuebingen_cohort_wo_celltype_correction//whole_cellcomposition_two_cohorts.pdf',width=16,height=9,l)




