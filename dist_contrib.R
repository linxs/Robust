dist.contrib<-function(mc=NULL,mre=NULL){
  
  eudist<-sqrt(((mc[,1,]-mre[,1,])^2)+((mc[,2,]-mre[,2,])^2))
  
  suma<-colSums(eudist)#sumatoria de distancias
  perc<-round((eudist/suma)*100,digits=4)#porcentual por punto
  
  results<-array(0,c(dim(mc)))
  results[,1,]<-round(eudist,digits=6)
  results[,2,]<-perc
  mvar<-mean(apply(eudist,2,var))
  
  print(paste("Residuals mean variance: ",round(mvar,digits = 6)),quote = F)
  
  print("Distance and % contribution to total (euclidean) distance between sides for each config:",quote = F)
  print(results)
  return(results)
}