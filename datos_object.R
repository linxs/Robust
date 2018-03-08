library(plyr)
v<-c(0,2,0,1,0,0,0,-1,0,-2,   -1.5,1, -1,1, -1.5,0, -1,0, -1.5,-1, -1,-1, 1,1, 1.5,1,1,0, 1.5,0, 1,-1, 1.5,-1.0   )


k<-2#dimensiones
sn<-10#numero de especimenes
ln<-17#numero de landmarks
sag<-5#n de sagitales
#lp<-c(2,6,9,17)
lp<-c(2,6,7,8,9)

pares<-as.matrix(read.csv("pairs.txt",header = F,sep=" "))



#B<-array(v,c(2,ln,sn))+rnorm((ln*2)*sn,sd=0.05)#descomentar para meter ruido de fondo

B<-array(v,c(k,ln,sn))
B[,lp,]<-B[,lp,]+rnorm(length(lp)*k*sn,sd=0.07)
plot(B[1,,],B[2,,])



#rotacion
alfa<-runif(sn,min = 0,max = pi/20)
AR<-sapply(alfa,simplify = "array",function(x){matrix(c(cos(x),-sin(x),sin(x),cos(x)),2,2)})
TAR<-array(0,c(k,ln,sn))
for(s in 1:sn){
  TAR[,,s]<-AR[,,s]%*%B[,,s]
}

#traslacion
TAR<-aaply(TAR,3,function(x){x+matrix(runif(2),k,ln)})
A<-aperm(TAR,c(3,2,1))
plot(A[,1,1],A[,2,1],asp = 1)

#-----------------------------------------

plot(A[,1,1],A[,2,1],asp = 1)
text(A[,1,1],A[,2,1],c(seq(1:(ln))),pos=2,offset=0.5,cex=0.55)

#pares<-matrix(c(6,13,7,12,8,15,9,14,10,17,11,16),6,2,byrow = T)


