#Esta funcion acepta un array A de dimensiones  pxkx2*n (n de landmarks, dimensiones, por ahora 2, y n de especimenes, cada uno con dos configuraciones)
#Las configuraciones deben estar ordenadas por individuo: lado I del ind 1, luego lado D del ind 1, etc...
#ctr puede tomar los valores "gmedian"(mediana espacial", "median" (mediana cac) y "mean" (media))

matching.symm<-function(A,ctr="gmedian",legend.loc="topleft"){
  
  library(MASS)
  library(Robgit)
  library(Gmedian)
 
  nl<-length(A[,1,1])#numero de landmarks
  n<-dim(A)[3]/2

  #--------------------------centrado----------------
  Mc<-center(A,cent=ctr)
 
  #---------- crea arrays por lado
  
  side<-c(1,2)
  Mcl<-Mc[,,side==1,drop=FALSE]
  Mcr<-Mc[,,side==2,drop=FALSE]
  
  dv<-Mcr-Mcl
  
  mod<-as.matrix(sqrt((dv[,1,]^2)+(dv[,2,]^2)))#norma
  N<-aperm(array(mod,c(dim(dv)[1],n,2)),c(1,3,2))#crea array con la norma repetida, para operar con el array de vectores (abajo)
  ndv<-dv/N #normaliza los vectores diferencia
  ndvp<-aperm(apply(ndv,c(1,3),function(x)if(x[2]<0) {x<-x*-1}else{x<-x}),c(2,1,3))#hace primera componente positiva
  par.mat<-as.matrix(acos(ndvp[,1,]))#calcula el angulo que corresponde al componente de cada vector sobre x

  
  dr<-apply(par.mat,2,median)#saca la dirección mediana (aún en radianes)
  vr<-cbind(cos(dr),sin(dr))#vectores unitarios definidos por sus proyecciones sobre x e y a partir de seno y coseno del angulo de menores proyecciones
  
  #calcula y aplica la matriz de Householder
  Ur<-array(0,dim(Mcl))
  for(i in 1:n){
    R<-diag(2) - 2*vr[i,]%*%t(vr[i,])
    Ur[,,i] <- Mcl[,,i]%*%R
  }
  
  S<-array(t(c(Mcr,Ur)),c(nl,2,(n*2)))
  print("performing final robust superimposition...",quote = F)
  capture.output(U<-robgit(S))#un ajuste robusto para eliminar pequeñas diferencias
  Ud<-U[,,1:n,drop=FALSE]
  Ui<-U[,,(n+1):(n*2),drop=FALSE]
  
  #saca contrib % de cada punto
  distances<-dist.contrib(mc=Ui,mre=Ud)
  
  #plots
  
  plot.result(mc=Ud, mre =Ui, nconf = n,legloc = legend.loc,object=F)
  return(list(distances,Ui-Ud))
}

