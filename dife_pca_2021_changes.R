
dife_pca <- function(input_df,filename='',colorgroups= NULL) {
  if (dim(input_df)[1] > 3) {
    require("PCAtools")

 
    #color coding
    if(is.null(colorgroups)){
      print("no groups setup, using standart colors")
      groups<-c(rep("1",ncol(input_df)))
    } else {
      print("groups setup, trying to use different colors per group")
      if(length(colorgroups)==ncol(input_df)){
        print("correct length, using the groups choosen")
        groups<-factor(colorgroups)
      } else {
        print("wrong number of groups, back to default")
        groups<-c(rep("1",ncol(input_df)))
      }
    }
    if(anyNA(input_df)){
      print("found NAs in dataset, we exchange them to zeros")
      input_df[is.na(input_df)] <- 0
    }
    
    metadata<-data.frame(row.names=colnames(input_df))
    metadata$group<-groups
    
    #PCA calculation
    result.pca <- pca(input_df,metadata = metadata)
    
    
  }
  
  return(result.pca)
}

