CONDSIMpreloop<-function(covParamObj,gridParamObj){

a<-Sys.time()

 ## source R functions
    setwd("/home/pwg/CONDSIM/CONDSIMalgorithm")
    source("makeFullPIxelLocationLists.R")
    source("PopulateCovarianceMatrices.R")
    source("calculateDistanceMatrices.R")
    source("getfootprintLUTs.R")
    source("getPixelSubset.R")
    source("getModifiedFootprintMatrices.R")
    source("getSubsetInterimMatrices.R")
    source("MVRNORM.R")

 ## source fortran functions
    system("R CMD SHLIB geographic.f")
    system("R CMD SHLIB euc_angle.f")
    system("R CMD SHLIB dist_to_aniso.f")
    system("R CMD SHLIB isotropic_cov_funs.f")
    system("R CMD SHLIB Stein_st_covariance_serial.f")
    dyn.load("geographic.so")
    dyn.load("euc_angle.so")
    dyn.load("dist_to_aniso.so")
    dyn.load("isotropic_cov_funs.so")
    dyn.load("Stein_st_covariance_serial.so")  
 
 ## invoke parameter file
    source("ParamFile_uncond.R")

 ## set up grid parameters    
    source("setupGridParams.R")

 ## calculate top edge position
    TOPEDGELAT<-YLLCORNER + (CELLSIZE*(Nrows))

 ## get lists of relative pixel locations (row,col,t) for all data and prediction pixels, and make data locatinos into standard table
    PixelLocationObj<-makeFullPIxelLocationLists()
    footprintXYT<-cbind(PixelLocationObj$xd,PixelLocationObj$yd,PixelLocationObj$td)
    yseq<-PixelLocationObj$yseq ## this is needed if we are thinning the footprint - defines which rows in each column are included
 
 ## use pixel locations to obtain absolute anisotropic distance matrices
    DistMatObj<-calculateDistanceMatrices(PixelLocationObj)

 ## populate covariance matrices
    CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
    
 ## define PtoP covariance matrix (this is constant - prediction set is always a single column - regardless of where we are predicting)   
    cPtoP<-CovMatObj$cPtoP

 ## loop through subset positions for full footprint and carry out matrix algebra for each. Assign output lists to first position of global lists
    SubsetInterimMatrices.Obj<-getSubsetInterimMatrices(CovMatObj$cDtoD,CovMatObj$cDtoP,cPtoP,returnL=T)
    PostMeanInterim.FULLlist<-SubsetInterimMatrices.Obj$PostMeanInterim.list
    PostVar.FULLlist<-SubsetInterimMatrices.Obj$PostVar.list

 ## create non-standard footprint LUTs (where full footprint is truncated due to being close to boundary in space or time):
    footprintLUTObj<-getfootprintLUTs(PixelLocationObj)
    footprintLUT<-footprintLUTObj$footprintLUT ## LUT for each month and col of which of these modified footprints is required
    modifiedFootprintTable<-footprintLUTObj$modifiedFootprintTable  ## table of unique modified footprints (as logical vectors of pixels in standard footprint to keep)
    NuniqueFPs<-footprintLUTObj$NuniqueFPs  ## how many unique footprints are we dealing with?

 ## for each unique modified footprint obtain the modified cDtoD and cDtoP covariance matrices
 ## but only if more than 1 modified footprint (first month/col is special case)
    if(NuniqueFPs>1){
      ModifiedFootprintMatricesObj<-getModifiedFootprintMatrices(footprintLUTObj,CovMatObj)
      cDtoDmodified.list<-ModifiedFootprintMatricesObj$cDtoDmodified.list
      cDtoPmodified.list<-ModifiedFootprintMatricesObj$cDtoPmodified.list

      ## for each unique modified footprint, loop through the row subsets and carry out matrix algebra
       
      # initialise lists that will hold a list of matrices for the different subsets for each unique modified footprint type
      PostMeanInterim.LISTS<-PostVar.LISTS<-vector("list")

      for(FP in 1:NuniqueFPs){
          SubsetInterimMatrices.Obj<-getSubsetInterimMatrices(cDtoDmodified.list[[FP]],cDtoPmodified.list[[FP]],cPtoP,returnL=T,FP)
          PostMeanInterim.LISTS[[FP]]<-SubsetInterimMatrices.Obj$PostMeanInterim.list
          PostVar.LISTS[[FP]]<-SubsetInterimMatrices.Obj$PostVar.list
      }
    }

 ## initialise list of length MonthDepth+1 to hold output matrices used in algorithm (and initialise first element as blank matrix to hold current month)
    OutMATlist<-vector("list",MonthDepth+1)
    OutMATlist<-lapply(OutMATlist,length) #this initialises each element of the list as a number (zero) rather than a NULL element - necessary later
    OutMATlist[[1]]<-matrix(nrow=Nrows,ncol=Ncols)
    
 ## bundle up objects to return in series of nested lists
    footprintObj<-list("footprintLUT"=footprintLUT,"footprintXYT"=footprintXYT,"modifiedFootprintTable"=modifiedFootprintTable,"yseq"=yseq)
    KrigMatObj<-list("cPtoP"=cPtoP,"PostMeanInterim.FULLlist"=PostMeanInterim.FULLlist,"PostVar.FULLlist"=PostVar.FULLlist,"PostMeanInterim.LISTS"=PostMeanInterim.LISTS,"PostVar.LISTS"=PostVar.LISTS)
    preLoopObj<-list("footprintObj"=footprintObj,"KrigMatObj"=KrigMatObj,"OutMATlist"=OutMATlist)    

 ## return this list
    return(preLoopObj)
}


CONDSIMmonthloop<-function(month,preLoopObj,OutMATlist){

 ## define inputs from bundled list
    footprintLUT<-preLoopObj$footprintObj$footprintLUT
    footprintXYT<-preLoopObj$footprintObj$footprintXYT
    modifiedFootprintTable<-preLoopObj$footprintObj$modifiedFootprintTable
    yseq<-preLoopObj$footprintObj$yseq

    cPtoP<-preLoopObj$KrigMatObj$cPtoP
    PostMeanInterim.FULLlist<-preLoopObj$KrigMatObj$PostMeanInterim.FULLlist
    PostVar.FULLlist<-preLoopObj$KrigMatObj$PostVar.FULLlist
    PostMeanInterim.LISTS<-preLoopObj$KrigMatObj$PostMeanInterim.LISTS
    PostVar.LISTS<-preLoopObj$KrigMatObj$PostVar.LISTS

 ## invoke parameter file
    source("ParamFile_uncond.R")

 ## set up grid parameters    
    source("setupGridParams.R")

 ## loop through each column of this month
    for(col in 1:Ncols){ #1:Ncols

if(VERBOSE==2)print(paste("month",month,"  col",col))
           
       ## initialise data for this col/month, and also vector to hold subset predictions
       Zdata<-NULL      

       ## calculate mean vector for this prediction column
       meanP<-rep(0,Nrows)
                 
       ## if first column and first month (or not using previous months), simulate directly from mvn using PtoP covariance matrix 
       if((month==1|MonthDepth==0) & col==1){
         OutMATlist[[1]][,1]<-MVRNORM(1,MU=meanP,COV=cPtoP)
       }

       ## otherwise identify footprint type, extract data, identify interim kriging matrices, and preform prediction
       if(!((month==1|MonthDepth==0) & col==1)){

         ## get footprint type for this col/month
         footprintType<-footprintLUT[month,col]

         # if the standard footprint:
         if(footprintType==0){ footprint<-footprintXYT}

         # otherwise make extractions from standard footprint using modifiedFootprintTable
         if(footprintType!=0){ footprint<-footprintXYT[modifiedFootprintTable[,footprintType],]}
 
         ## extract data for this col/month given footprint type
         for(i in 0:min(c((month-1),MonthDepth))){
            dataCols<-unique(col+footprint[footprint[,3]==-i,1])
            if(class(THIN)!="NULL"){
              Zdata<-c(Zdata,as.vector(OutMATlist[[i+1]][yseq,dataCols]))
            }

            if(class(THIN)=="NULL"){
              Zdata<-c(Zdata,as.vector(OutMATlist[[i+1]][,dataCols]))
            }
         } 

         ## identify correct list of interim kriging matrices for this footprint type
         # if the standard footprint:
         if(footprintType==0){
           PostMeanInterim.list<-PostMeanInterim.FULLlist
           PostVar.list<-PostVar.FULLlist
         }     

         # if not the standard footprint:
         if(footprintType!=0){
           PostMeanInterim.list<-PostMeanInterim.LISTS[[footprintType]]
           PostVar.list<-PostVar.LISTS[[footprintType]]
         } 

         ## define mean vector for this data footprint
         meanD<-rep(0,length(Zdata))

         ## loop through the NSUBSETS subsets for this prediction column
         for(i in 1:NSUBSETS){

            ## which rows in the prediction column are we predicting?
            predIndex<-((i*SUBSETCOL)-(SUBSETCOL-1)):(i*SUBSETCOL)

            ## redfine mean vector under predictions for this subset
            meanP.subset<-meanP[predIndex]

            ## if first subset, no aditional data from this column is available
            if(i==1){

              ## perform remaining matrix algebra to obtain prediction vector for this subset of the prediction column
              PostMean.subset<-meanP.subset + PostMeanInterim.list[[i]]%*% (Zdata-meanD) 
            }

            ## if any other subset, we need to include data from SUBSETCOL already predicted pixels in the prediction column
            if(i>1){

              ## which rows in the prediction column have already been predicted and we are using as data?
              dataIndex<-predIndex-SUBSETCOL

              ## redefine Zdata and mean under data
              Zdata.subset<-c(Zdata,heldPredictedSubset)
              rm(heldPredictedSubset)

              meanD.subset<-c(meanD,meanP[dataIndex])

              ## perform remaining matrix algebra to obtain prediction vector for this subset of the prediction column
              PostMean.subset<-meanP.subset + PostMeanInterim.list[[i]]%*% (Zdata.subset-meanD.subset) 
            }

            ## now draw a realisation of this mvn vector and add this as a new subset of data for this column of the output matrix 
            temp<-try(heldPredictedSubset<-MVRNORM(1,MU=PostMean.subset,L=PostVar.list[[i]]))
            if(class(temp)=="try-error"){print(paste("!!!!!!!!!!!!!WARNING: MVRNORM FAILED FOR MONTH",month,", col",col))}
            OutMATlist[[1]][predIndex,col]<-heldPredictedSubset
         }
       }
    } 

 ## define grid to return for this month
    MonthGrid<-OutMATlist[[1]]
    
 ## this current month now becomes month t-1, other MonthDepth months all shuffle along
  # loop through matrices in memory and shuffle them all back in time by one (if we are storing any previous months)
    if(MonthDepth>0){ 
      for(i in MonthDepth:1){ 
         OutMATlist[[i+1]]<-OutMATlist[[i]]
      }
    }

 ## initialise a new blank matrix for the next current month
    OutMATlist[[1]]<-matrix(nrow=Nrows,ncol=Ncols)
 
 ## return grid for this month and OutMATlist object to pass on to subsequent months
    monthObject<-list("MonthGrid"=MonthGrid,"OutMATlist"=OutMATlist)
    return(monthObject)
}    


