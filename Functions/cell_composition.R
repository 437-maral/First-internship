cell_composition=function(data,cell,phenotypes,names3,names5,names6){
  data=select(data,cell)
  data=min_max_norm(data)
  
  rownames(data)=phenotypes$Sample_Name
  rownames(data)[which(rownames(data) %in% names3)]<-paste0(rownames(data)[which(rownames(data) %in% names3)],"_3")
  rownames(data)[which(rownames(data) %in% names5)]<-paste0(rownames(data)[which(rownames(data) %in% names5)],"_5")
  rownames(data)[which(rownames(data) %in% names6)]<-paste0(rownames(data)[which(rownames(data) %in% names6)],"_6")
  rownames(data)
  
  data_t=t(data)
  return(data_t)
}
