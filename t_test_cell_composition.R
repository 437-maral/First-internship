####doing t_test
##data=cell_composition

cell_com=fread('L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/cellcomposition.tsv')




cell_com=select(cell_com,-Name)
cell_com$Group=as.character(cell_com$Group)


###compare between clusters

plot=list()
for(i in 2:ncol(cell_com)-1){
  plot[[i]]=
    ggplot(cell_com,aes_string(x = 'Group',y = names(cell_com)[i],fill='Group')) +
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
    ggtitle(names(cell_com)[i])
}


ggsave(paste0("L:/Bioinformatics_Unit/T2D_Subcluster_Studie/potsdam_cohort_wo_celltype_correction/cell_composition", 'cell_composotion_potsdam', ".pdf"),width=16,height=9,  gridExtra::marrangeGrob(grobs = plot, nrow = 2, ncol = 3))

	
