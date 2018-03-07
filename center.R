center<-function(A,cent){
  library(Gmedian)
    ni<-dim(A)[3]
    switch(cent,
         mean={Mc<-apply(A,c(2,3), function(x) x-mean(x))},
         median={Mc<-apply(A,c(2,3), function(x) x-median(x))},
         gmedian={Mc<-array(0,c(dim(A)))
         gmed<-t(apply(A,c(3), function(x) Weiszfeld(x)))
         for(i in 1:ni){
           Mc[,,i]<-t(t(A[,,i])-c(gmed[[i]][[1]]))
         }    
         }
  )
  return(Mc)
}
