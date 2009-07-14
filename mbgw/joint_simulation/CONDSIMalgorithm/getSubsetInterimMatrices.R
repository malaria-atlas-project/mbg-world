## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################

getSubsetInterimMatrices<-function(cDtoD,cDtoP,cPtoP,returnL,NuniqueFPs,counter=0){

# cDtoD<-CovMatObj$cDtoD;cDtoP<-CovMatObj$cDtoP
# cDtoD<-ModifiedFootprintMatricesObj$cDtoDmodified.list[[2]];cDtoP<-ModifiedFootprintMatricesObj$cDtoPmodified.list[[2]];counter<-2

 ## if we are dealing with the first unique modified footprint, this is a special case (first col/month - will be simulated directly)
 ## we do not need to define subset matrices so return null for this list location
    if(counter==1){
      if(VERBOSE>=1) print(paste("Modified FP",counter,"of",NuniqueFPs,": no matrices to solve"))
      return(list("PostMeanInterim.list"=NULL,"PostVar.list"=NULL)) 
    }

 ## initialise lists in which matrices output will be stored
    PostMeanInterim.list<-vector("list")
    PostVar.list<-vector("list")

 ## define number of data and predictions in this full footprint
    Ndata<-ncol(cDtoD)
    Npred<-ncol(cPtoP)

 ## define combined matrices from which all subst matrices can be extracted
    cDtoDcombined<-rbind(cbind(cDtoD,cDtoP),cbind(t(cDtoP),cPtoP))
    cDtoPcombined<-rbind(cDtoP,cPtoP)

 ## loop through each of the NSUBSETS subset positions and define appropriate matrices
    for(i in 1:NSUBSETS) {

    if(VERBOSE>=1){
      if(counter==0)print(paste("Full FP : solving for subset",i,"of",NSUBSETS))
      if(counter!=0)print(paste("Modified FP",counter,"of",NuniqueFPs,": solving for subset",i,"of",NSUBSETS))
    }        
       ## which rows in the prediction column are we predicting?
       predIndex<-((i*SUBSETCOL)-(SUBSETCOL-1)):(i*SUBSETCOL)

       ## modify matrices for first subset of column (differs because no pixels in pred column are data)
       if(i==1){

         ## define subseted cDtoD matrix as original cDtoD matrix (no data from prediction column)
         cDtoD.temp<-cDtoD

         ## define subseted cDtoP matrix: remove those pixels from prediction column not being pedicted this subset
         cDtoP.temp<-cDtoP[,predIndex]  

         ## define subseted cPtoP matrix to include only those pixels in prediction row that sre being predicted this subset
         cPtoP.temp<-cPtoP[predIndex,predIndex] 
       }

       ## modify matrices for all other subsets of column (SUBSETCOL pixels preceding subset in pred column are data)
       if(i>1){
         ## and which rows in the prediction column are we using as data, having already been predicted?
         dataIndex<-predIndex-SUBSETCOL

         ## define subseted cDtoD matrix to include aditional data from prediction row
         cDtoD.temp<-cDtoDcombined[c(1:Ndata,dataIndex+Ndata),c(1:Ndata,dataIndex+Ndata)]  

         ## define subseted cDtoP matrix: some pixels from prediction column are data, others are predictions
         cDtoP.temp<-cDtoPcombined[c(1:Ndata,dataIndex+Ndata),predIndex]  

         ## define subseted cPtoP matrix to include only those pixels in prediction row that are being predicted this subset
         cPtoP.temp<-cPtoP[predIndex,predIndex] 
       }

       ## carry out matrix algebra on the modified matrices for this subset
          cDtoD.inv<-solve(cDtoD.temp)
          PostMeanInterim<-t(cDtoP.temp) %*% cDtoD.inv
          PostVar<-cPtoP.temp - t(cDtoP.temp) %*% cDtoD.inv %*% cDtoP.temp

       ## if we are returning the choleski decomposition of the kriging variance matrix, rather than the original matrix, do this now:
          if(returnL){
       
            ichol_full.list<-.Fortran("ichol_full",              
                       c=as.double(PostVar),              
                       n=as.integer(Npred),              
                       sig=as.double(rep(0,Npred^2)),              
                       m=as.integer(0),              
                       p=as.integer(rep(0,Npred)))              
            U<-matrix(ichol_full.list$sig,nrow=Npred,ncol=Npred)              
            n<-ichol_full.list$m
            pivot<-ichol_full.list$p

#            U<-chol(PostVar, pivot = TRUE)
#            pivot <- attr(U, "pivot")
#            n <-attr(U,"rank")
            oo <- order(pivot)
            L<-t(U[1:n,oo])
            PostVar<-L 
          }       
          print(paste("range of L from getSubsetInterimMatrices:",min(L),"to",max(L)))
          if(any(is.nan(L))) print("ERROR!!! in getSubsetInterimMatrices, found nan values in L") 
          if(any(is.na(L))) print("ERROR!!! in getSubsetInterimMatrices, found na values in L") 
          if(any(is.infinite(L))) print("ERROR!!! in getSubsetInterimMatrices, found Inf values in L")       
      

       ## add outputs to relevant lists
          PostMeanInterim.list[[i]]<-PostMeanInterim 
          PostVar.list[[i]]<-PostVar
          rm(PostMeanInterim)           
          rm(PostVar)
   }

 ## return matrices lists for this footprint 
    return(list("PostMeanInterim.list"=PostMeanInterim.list,"PostVar.list"=PostVar.list)) 
}
