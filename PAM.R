
pam_function <- function(input,title="PAM",class="pheno"){

pa <-pam(t(input), 3, diss =  FALSE,
                    metric = "manhattan", 
                    medoids =  "random",
                   # nstart =  10 ,
                    nstart =  100 ,
                    stand = FALSE, cluster.only = FALSE,
                    do.swap = TRUE,
                    keep.diss = FALSE,
                    keep.data = FALSE,
                    variant = "original",
                    trace.lev = 0)

#to create a confusion matrix for the 
clustering <- as.data.frame(pa[["clustering"]])
clustering <- rownames_to_column(clustering)
#get only the number of the cluster
clustering$rowname <- gsub(".+_","", clustering$rowname) %>% as.factor()
colnames(clustering) <-c("target","prediciton")
#important to change 3 first because it is in both target and prediciton
#pam would always say 1,2,3 as cluster, but we have 3,5,6 therefore we change here

if(class=="sofya"){#phenos in 8_pam, all others in pam function which are not MEP samples
   clustering$prediciton <- gsub(2,6,clustering$prediciton) 
   clustering$prediciton <- gsub(3,5,clustering$prediciton)
   clustering$prediciton <- gsub(1,3,clustering$prediciton)
} else if(class=="markus") { #others in pam_function of sofya, she used that for the MEP control data
   clustering$prediciton <- gsub(3,6,clustering$prediciton) 
   clustering$prediciton <- gsub(2,5,clustering$prediciton) 
   clustering$prediciton <- gsub(1,3,clustering$prediciton) 
} else if(class=="unknown1") {
  clustering$prediciton <- gsub(1,5,clustering$prediciton) 
  clustering$prediciton <- gsub(3,3,clustering$prediciton)
  clustering$prediciton <- gsub(2,6,clustering$prediciton)  
} else if(class=="unknown2") {
  clustering$prediciton <- gsub(1,6,clustering$prediciton) 
  clustering$prediciton <- gsub(3,5,clustering$prediciton)
  clustering$prediciton <- gsub(2,3,clustering$prediciton)  
} else if(class=="unknown3") {
  clustering$prediciton <- gsub(1,5,clustering$prediciton) 
  clustering$prediciton <- gsub(3,6,clustering$prediciton)
  clustering$prediciton <- gsub(2,3,clustering$prediciton)  
} else if(class=="unknown4") {
  clustering$prediciton <- gsub(1,6,clustering$prediciton) 
  clustering$prediciton <- gsub(2,5,clustering$prediciton)
  clustering$prediciton <- gsub(3,3,clustering$prediciton)  
}


 
 
 

pa$data <- t(input)
#vizualize
l <- fviz_cluster(pa,ellipse.type = "t", # Concentration ellipse
                             repel = TRUE, # Avoid label overplotting (slow)
                             main = title,
                             palette = "Set2",
                             pointsize = 8,
                             shape = 16,
                             geom = "point",
                             legend.title = "Predicted Cluster",
                             show.legend.text = FALSE,
                             ggtheme =  theme_classic(base_size = 20))+
  geom_text(aes(label = clustering$target),size=7,color="black")+
  #geom_point(aes(shape = clustering$target))
  scale_fill_discrete(labels=c("3","6","5")) +
  scale_color_discrete(labels = c("3","6","5"))+
  guides(col = FALSE)
print(l)

conf_mat <- confusion_matrix(targets = clustering$target,
                             predictions = clustering$prediciton)

print(plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]],
                            font_counts = font(size=16),
                            add_row_percentages = FALSE,
                            add_col_percentages = FALSE,
                            arrow_size = 0.1,
                            palette = "Greens",
                            rotate_y_text = F,
                            add_normalized = FALSE)+theme_minimal(base_size = 28))

d <- table(clustering)
acc <- sum(diag(d))/sum(d) #overall accuracy
inccorects <- 1-sum(diag(d))/sum(d) #incorrect classification 

print(paste0(
  "Samples = ", ncol(input), "  Features = ", nrow(input), 
  "  Silhouette = ",round(pa$silinfo$avg.width,2), "  Accuracy = ",round(acc*100,2))
)
}
