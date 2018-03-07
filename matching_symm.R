#Esta funcion acepta un array A de dimensiones  pxkx2*n (n de landmarks, dimensiones, por ahora 2, y n de especimenes, cada uno con dos configuraciones)
#Las configuraciones deben estar ordenadas por individuo: lado I del ind 1, luego lado D del ind 1, etc...
#Lee del dir de trabajo un archivo (por defecto llamado "side.txt") que es una lista del lado de cada config de A (DEBE usar L y R para izquierda y derecha)
#ctr puede tomar los valores "gmedian"(mediana espacial", "median" (mediana cac) y "mean" (media))

matching.symm<-function(A,ctr="gmedian",side.file="side.txt",legend.loc="topleft"){
  
  library(MASS)
  library(Robgit)
  library(Gmedian)
  
  side<-read.table(side.file,sep=" ",header=F)
  nd<-length(A[1,,1])#numero de dimensiones
  nl<-length(A[,1,1])#numero de landmarks
  n<-dim(A)[3]/2

  #--------------------------centrado----------------
  Mc<-center(A,cent=ctr)
 
  #---------- crea arrays por lado
  
  Mcl<-Mc[,,side=="L"]
  Mcr<-Mc[,,side=="R"]
  
  dv<-Mcl-Mcr
  
  ang.par<-atan(dv[,2,]/dv[,1,])#angulo en radianes
  par.mat<-apply(ang.par,c(1,2),function(x)if(x<0) {pi-abs(x)}else{x})#esto hace que los ángulos sean todos positivos entre 0 y pi
  
  dr<-apply(par.mat,2,median)#saca la dirección mediana (aún en radianes)
  vr<-cbind(cos(dr),sin(dr))#vectores unitarios definidos por sus proyecciones sobre x e y a partir de seno y coseno del angulo de menores proyecciones
  
  #calcula y aplica la matriz de Householder
  Ur<-array(0,dim(Mcl))
  for(i in 1:n){
    R<-diag(2) - 2*vr[i,]%*%t(vr[i,])
    Ur[,,i] <- Mcl[,,i]%*%R
  }
  
  S<-array(t(c(Mcr,Ur)),c(nl,2,(n*2)))
  print("please wait...",quote = F)
  sink("/dev/null")
  U<-robgit(S)#un ajuste robusto para eliminar pequeñas diferencias
  sink()
  Ui<-U[,,1:n]
  Ud<-U[,,(n+1):(n*2)]
  
  #saca contrib % de cada punto
  distances<-dist.contrib(mc=Ui,mre=Ud)
  
  #plots
  plot.result(mc=Ud, mre =Ui, nconf = n,legloc = legend.loc,object=F)
  return(list(distances,Ui-Ud))
}

