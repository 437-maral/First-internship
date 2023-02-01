####doing t_test
##data=cell_composition


setwd('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction')
cell_com=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/cellcomposition.tsv')




cell_com=select(cell_com,-Name)
cell_com$Group=as.character(cell_com$Group)

###my comparison

my_comparisons <- list(c("3", "5"),c("3", "6"),c("5", "6"))

##I need to speerate coulmn
cell_com3=cell_com[cell_com$Group==3,]
cell_com5=cell_com[cell_com$Group==5,]
cell_com6=cell_com[cell_com$Group==6,]


cell_com3_5=rbind(cell_com3,cell_com5)
cell_com3_6=rbind(cell_com3,cell_com6)
cell_com5_6=rbind(cell_com5,cell_com6)

###compare between clusters
###comparison between 3_5


plot3_5=list()
for(i in 2:ncol(cell_com3_5)-1){
  plot3_5[[i]]=
    ggplot(cell_com3_5,aes_string(x = 'Group',y = names(cell_com3_5)[i],fill='Group')) +
    geom_boxplot(show.legend = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
    stat_compare_means(label.x.npc="center")+
    theme_classic()+
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          #axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ggtitle(names(cell_com3_5)[i])+
    scale_fill_brewer(palette="Dark2")
}


###comparison 3_6


plot3_6=list()
for(i in 2:ncol(cell_com3_6)-1){
  plot3_6[[i]]=
    ggplot(cell_com3_6,aes_string(x = 'Group',y = names(cell_com3_6)[i],fill='Group')) +
    geom_boxplot(show.legend = FALSE)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
    stat_compare_means(label.x.npc="center")+
    theme_classic()+
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          #axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ggtitle(names(cell_com3_6)[i])
}

plot5_6=list()
for(i in 2:ncol(cell_com5_6)-1){
  plot5_6[[i]]=
    ggplot(cell_com5_6,aes_string(x = 'Group',y = names(cell_com5_6)[i],fill='Group')) +
    geom_boxplot(show.legend = FALSE)+
    scale_fill_manual(values=c("#9ac0cd","#9bcd9a"))+
    geom_point(position=position_jitterdodge(jitter.width =0.1),show.legend = FALSE)+
    stat_compare_means(label.x.npc="center")+
    theme_classic()+
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          #axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ggtitle(names(cell_com5_6)[i])
}

pdf('whole_cellcomposition.pdf',width=16,height=9)
gridExtra::marrangeGrob(grobs = plot3_5, nrow = 2, ncol = 3,top='cluster3_5')
gridExtra::marrangeGrob(grobs = plot3_6, nrow = 2, ncol = 3,top='cluster3_6')
gridExtra::marrangeGrob(grobs = plot5_6, nrow = 2, ncol = 3,top='cluster5_6')
dev.off()

