#Esta funcion acepta un array A de dimensiones  pxkxn (n de landmarks, dimensiones, por ahora 2, y n de especimenes)
#Los landmarks de cada configuración deben estar en el siguiente orden: sagitales, izquierdos y derechos (o der e izq, eso es indiferente)
#Lee del dir de trabajo un archivo (por defecto llamado "pairs.txt") que es una lista de los pares de landmarks
#ctr puede tomar los valores "gmedian"(mediana espacial", "median" (mediana cac) y "mean" (media))
#opciones referencias: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"

object.symm<-function(A,ctr="gmedian",prs.file="pairs.txt",proj.met="msum",legend.loc="topleft"){
  
  library(Gmedian)
  
  pares<-read.table(prs.file,sep=" ",header=FALSE)
  ns<-length(A[,1,1])-(length(pares[,1])*2)#numero de puntos sagitales
  np<-length(A[,1,1])-ns #numero de puntos pareados
  n<-dim(A)[3]
  
  #--------------------------centrado----------------
  
  Ac<-center(A,cent=ctr)
  
  #lee la base y crea matrices para sagitales y para cada lado, izquierdo y derecho (en el txt están ordenados así)
  m0<-Ac[1:ns,,,drop=FALSE]
  m1<-Ac[(ns+1):((ns+1)+(np/2)-1),,,drop=FALSE]
  m2<-Ac[((ns+1)+(np/2)):((ns)+(np)),,,drop=FALSE]
  
  #puntos sagitales
  pares.sag<-t(combn(1:ns,2))# aca hace la combinatoria del n de orden de los sagitales tomados de a 2
  
  
  vs<-m0[pares.sag[,1],,,drop=FALSE]-m0[pares.sag[,2],,,drop=FALSE]#obtiene los vectores de direcciones sagitales para los pares.sag
  
  mod<-as.matrix(sqrt((vs[,1,]^2)+(vs[,2,]^2)))#norma
  N<-aperm(array(mod,c(dim(pares.sag)[1],n,2)),c(1,3,2))#crea array con la norma repetida, para operar con el array de vectores (abajo)
  vsn<-vs/N#vectores sagitales normalizados

  
  #puntos pareados
  vp<-m1[pares[,1]-ns,,,drop=FALSE]-m2[pares[,2]-ns-(np/2),,,drop=FALSE]
  
  # proyecta y ordena las direcciones sagitales por proyeccion

  
  P<-array(0,dim = c(dim(vs)[1],(np/2),n))#crea un array para guardar las proyecciones
  
  #Abajo, proyecta sacando producto punto entra la matriz de vectores sagitales normalizados y
  #la de vectores pareados para cada "feta" del array correspondiente (cada configuracion) 
 
  for(i in 1:n){
    P[,,i]<-(vsn[,,i]%*%t(vp[,,i]))
  }
  
  P<-abs(P)
  
  p.sum<-apply(P,c(1,3), sum)#suma de proyecciones. Cada columna es un especimen. Adentro, cada fila es la suma de proyecciones de los pareados
  p.median<-apply(P,c(1,3),median)#mediana de proyecciones
  
  # indices del vector sagital con suma o mediana de proyecciones minima para cada configuracion
  switch(proj.met,
    msum={i.min<-apply(p.sum,c(2),which.min)},
    mmedian={i.min<-apply(p.median,c(2),which.min)}
  )
  
  vr<-NULL #crea matriz para vectores sagitales unitarios para reflexion (que luego se completa abajo)
  
  #Selecciona los vectores normalizados que tuvieron proyeccion minima y los copia a la matriz vr
  for(j in 1:n){
    vaux<-vsn[i.min[j],,j]
    vr<-rbind(vr,vaux)
  }
  #--------------------------reflexion----------------
  
  #transforma al vector en ortogonal al eje de reflexion (es decir al vector que unia los 2 puntos)
  e<-matrix(0,n,2)
  e[,2]<-vr[,1]
  e[,1]<- -(vr[,2])
  
  #calcula y aplica la matriz de Householder
  Ur<-array(0,dim(A))
  for(i in 1:n){
    R<-diag(2) - 2*e[i,]%*%t(e[i,])
    Ur[,,i] <- Ac[,,i]%*%R
  }
  
  #reetiquetado
  Ure<-Ur #hace un duplicado de la matriz rotada, y el for de abajo cambia los izquierdos por los derechos y viceversa 
  for (j in 1:(np/2)) {
    Ure[pares[j, 1],,] <- Ur[pares[j,2],,]
    Ure[pares[j, 2],,] <- Ur[pares[j,1],,]
  }
  
  T<-(Ac-Ure)/2 #saca la diferencia entre original y reflejada y lo divide por 2 (lo que queda,T, es el residuo entre original y simétrica)
  
  #saca contrib % de cada punto 
  distances<-dist.contrib(mc=Ac,mre=Ure)
  #plots
  plot.result(mc=Ac,mre=Ure,mt=T,nconf=n,object=TRUE,legloc=legend.loc)
  #devuelve T
  return(list(distances,T))
  
}
