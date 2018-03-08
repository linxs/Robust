
plot.result<-function(mc=NULL,mre=NULL,mt=NULL,nconf,object=TRUE,legloc="topleft"){
  conf<-1:nconf
  if(object[1]){
    for(c in conf){
      plot(mc[,,c],asp=1,xlab = "",ylab = "",main = paste("Config",c))#original
      points(0,0,pch=3)#centro (0,0)
      points(mre[,,c],pch=20)#reflejada
      points((mt[,,c]+mre[,,c]),col="red",pch=20,cex=0.8)#simétrica (residuos + reflejada)
      text(mc[,1,c],mc[,2,c],c(seq(1:(dim(mc)[1]))),pos=2,offset=0.5,cex=0.55)#numero de landmark
      legend(x=legloc,c("Original","Reflected","Symmetric"),pch=c(1,20,20) ,col = c("black","black","red"), cex=0.7)
    }
    
  }else{
    
    for(c in conf){
      plot(mc[,,c],asp=1,xlab = "",ylab = "",main = paste("Config",c))#original
      points(0,0,pch=3)#centro (0,0)
      points(mre[,,c],col="red",pch=20,cex=0.8)#reflejada
      text(mc[,1,c],mc[,2,c],c(seq(1:dim(mc)[1])),pos=2,offset=0.5,cex=0.55)
      legend(x=legloc,c("Left","Right"),pch=c(20,1) ,col = c("red","black"), cex=0.7)
    }
  }
  
}
  
