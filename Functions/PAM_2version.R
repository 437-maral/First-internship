####just pam function

library(data.table)
library(cluster)
library(tibble)
library(tidyr)
library(gridExtra)
library(factoextra)
library(ggrepel)
library(cluster)
library(dplyr)
library(stringr)
library(openxlsx)
library(cvms)


###Functions

Pam_function<-function(dataset,name){
  set.seed(666)
  pa <-pam(t(dataset), 3, diss =  FALSE,
           metric = "manhattan", 
           medoids =  "random",
           nstart =  10 ,
           stand = FALSE, cluster.only = FALSE,
           do.swap = TRUE,
           keep.diss = FALSE,
           keep.data = FALSE,
           variant = "faster",
           trace.lev = 0)
  
  pa$data <- t(dataset)
  clustering <- as.data.frame(pa[["clustering"]])
  clustering <- rownames_to_column(clustering)
  #get only the number of the cluster
  clustering$rowname <- gsub(".+_","", clustering$rowname) %>% as.factor()
  colnames(clustering) <-c("target","prediciton")
  plots<-list()
  plots[[1]]<-confmatrix_plotter(clustering,class="sofya",name)
  plots[[2]]<-confmatrix_plotter(clustering,class="markus",name)
  plots[[3]]<-confmatrix_plotter(clustering,class="unknown1",name)
  plots[[4]]<-confmatrix_plotter(clustering,class="unknown2",name)
  plots[[5]]<-confmatrix_plotter(clustering,class="unknown3",name)
  plots[[6]]<-confmatrix_plotter(clustering,class="unknown4",name)
  
  ggsave(filename=paste0(main_folder,name,"_PAM_conf.pdf"),width = 12, height = 12,gridExtra::marrangeGrob(grobs = plots, nrow = 2, ncol = 3))
  
  
  cat(paste0("PAM for: ",name," done \n"))
}

confmatrix_plotter<-function(clustering_data,class,name){
  #important to change 3 first 
  #pam would always say 1,2,3 as cluster, but we have 3,5,6 therefore we change here
  if(class=="sofya"){
    clustering_data$prediciton <- gsub(2,6,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(3,5,clustering_data$prediciton)
    clustering_data$prediciton <- gsub(1,3,clustering_data$prediciton)
  } else if(class=="markus") { 
    clustering_data$prediciton <- gsub(3,6,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(2,5,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(1,3,clustering_data$prediciton) 
  } else if(class=="unknown1") {
    clustering_data$prediciton <- gsub(1,5,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(3,3,clustering_data$prediciton)
    clustering_data$prediciton <- gsub(2,6,clustering_data$prediciton)  
  } else if(class=="unknown2") {
    clustering_data$prediciton <- gsub(1,6,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(3,5,clustering_data$prediciton)
    clustering_data$prediciton <- gsub(2,3,clustering_data$prediciton)  
  } else if(class=="unknown3") {
    clustering_data$prediciton <- gsub(1,5,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(3,6,clustering_data$prediciton)
    clustering_data$prediciton <- gsub(2,3,clustering_data$prediciton)  
  } else if(class=="unknown4") {
    clustering_data$prediciton <- gsub(1,6,clustering_data$prediciton) 
    clustering_data$prediciton <- gsub(2,5,clustering_data$prediciton)
    clustering_data$prediciton <- gsub(3,3,clustering_data$prediciton)  
  }
  d <- table(clustering_data)
  acc <- sum(diag(d))/sum(d) #overall accuracy
  
  conf_mat <- confusion_matrix(targets = clustering_data$target,
                               predictions = clustering_data$prediciton)
  
  
  g<-plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]],
                           font_counts = font(size=16),
                           add_row_percentages = FALSE,
                           add_col_percentages = FALSE,
                           arrow_size = 0.1,
                           palette = "Reds",
                           rotate_y_text = F,
                           add_normalized = FALSE)+
    theme_minimal(base_size = 26)+
    theme(plot.title = element_text(size=16))+
    ggtitle(paste0(class, " acc:",round(acc*100,2),"%"))
  
  
  return(g)
  
  
}


pheno_preper <- function (choosen_pheno,names3,names5,names6){
  #we only can take phenotypes where we have data for ALL samples, kick the rest
  choosen_pheno <- choosen_pheno[rowSums(is.na(choosen_pheno)) == 0,]
  rownames(choosen_pheno) <- choosen_pheno$sample_name
  #normalise with minmax normalisation
  
  cluster_choosen_pheno_norm <-  as.data.frame(lapply(choosen_pheno[2:ncol(choosen_pheno)], min_max_norm))
  cluster_choosen_pheno_norm_t <- t(cluster_choosen_pheno_norm)
  colnames(cluster_choosen_pheno_norm_t) <- choosen_pheno$sample_name
  colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names3] <- paste0(colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names3] ,"_3")
  colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names5] <- paste0(colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names5] ,"_5")
  colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names6] <- paste0(colnames(cluster_choosen_pheno_norm_t)[colnames(cluster_choosen_pheno_norm_t)%in% names6] ,"_6")
  
  return(cluster_choosen_pheno_norm_t)
}