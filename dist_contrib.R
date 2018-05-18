dist.contrib<-function(mc=NULL,mre=NULL){
  
  eudist<-as.matrix(sqrt(((mc[,1,]-mre[,1,])^2)+((mc[,2,]-mre[,2,])^2)))
  
  suma<-colSums(eudist)#sumatoria de distancias
  perc<-round((eudist/suma)*100,digits=4)#porcentual por punto
  
  results<-array(0,c(dim(mc)))
  results[,1,]<-round(eudist,digits=6)
  results[,2,]<-perc
  
  print("Distance and % contribution to total (euclidean) distance between sides for each config:",quote = F)
  print(results)
  return(results)
}