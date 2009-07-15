MVRNORM<-function(Ndraws,MU,COV=c(),L=c()){

#Ndraws<-1;MU=PostMean.local;COV=PostVar.local

 # if we are passing the covriance matrix, need to define L
   if(class(L)=="NULL"){
   
     print("passed a covariance matrix, so performing cholesky decomp")


     ichol_full.list<-.Fortran("ichol_full",              
           c=as.double(COV),              
           n=as.integer(nrow(COV)),              
           sig=as.double(rep(0,nrow(COV)^2)),              
           m=as.integer(0),              
           p=as.integer(rep(0,nrow(COV))))              
    U<-matrix(ichol_full.list$sig,nrow=nrow(COV),ncol=nrow(COV))              
    n<-ichol_full.list$m
    pivot<-ichol_full.list$p

#            U<-chol(PostVar, pivot = TRUE)
#            pivot <- attr(U, "pivot")
#            n <-attr(U,"rank")
    oo <- order(pivot)
    L<-t(U[1:n,oo])
    rm(U)
    print(paste("range of L from MVRNORM:",min(L),"to",max(L)))

     # define choleski decomposition of COV
#     U<-chol(COV, pivot = TRUE)
#     pivot <- attr(U, "pivot")
#     n <-attr(U,"rank")
#     oo <- order(pivot)
#     L<-t(U[1:n,oo])

    #print(paste("n=",n))
    #print(paste("nrow of L=",nrow(L)))
    #print(paste("ncol of L=",ncol(L)))
   
   }
   
   
    
   #check for any NAs, nans, or Infs
   if(any(is.nan(L))) print("ERROR!!! in MVRNORM, found nan values in L") 
   if(any(is.na(L))) print("ERROR!!! in MVRNORM, found na values in L") 
   if(any(is.infinite(L))) print("ERROR!!! in MVRNORM, found Inf values in L")       
   
   n<-ncol(L)

 # take NDraws samples from the multivariate normal distribution of mean MU and covariance COV
   samples <- as.vector(MU) + (L %*% matrix(rnorm(n*Ndraws),nrow=n,ncol=Ndraws))
   rm(L)
   #print(paste("range of samples from MVRNORM:",min(samples),"to",max(samples)))

   return(t(samples))
}
