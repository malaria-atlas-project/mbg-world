## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################


CONDSIMpreloop<-function(covParamObj,gridParamObj,monthParamObj,startRel,endRel, paramfileINDEX){

#CONDSIMpreloop<-function(covParamObj,gridParamObj,monthParamObj){

    covParamObj <<- as.list(covParamObj)
    gridParamObj <<- as.list(gridParamObj)
    monthParamObj <<- as.list(monthParamObj)
    
       
 ## source R functions
#   setwd("/home/pwg/CONDSIM/CONDSIMalgorithm")
    source("makeFullPIxelLocationLists.R")
    source("PopulateCovarianceMatrices.R")
    source("calculateDistanceMatrices.R")
    source("getfootprintLUTs.R")
    source("getPixelSubset.R")
    source("getModifiedFootprintMatrices.R")
    source("getSubsetInterimMatrices.R")
    source("describeList.R")    
    source("MVRNORM.R")

 ## source fortran functions
    temp<-system("R CMD SHLIB geographic.f",intern=T)
    temp<-system("R CMD SHLIB euc_angle.f",intern=T)
    temp<-system("R CMD SHLIB dist_to_aniso.f",intern=T)
    temp<-system("R CMD SHLIB isotropic_cov_funs.f",intern=T)
    temp<-system("R CMD SHLIB Stein_st_covariance_serial.f",intern=T)
    dyn.load("geographic.so")
    dyn.load("euc_angle.so")
    dyn.load("dist_to_aniso.so")
    dyn.load("isotropic_cov_funs.so")
    dyn.load("Stein_st_covariance_serial.so")  
  
 ## invoke parameter file
    source(paste("ParamFile_uncond_",paramfileINDEX,".R",sep=""))
#print(paste("ColDepth=",ColDepth))
#print(paste("MonthDepth=",MonthDepth))
#print(paste("THINX=",THINX))
#print(paste("THINY=",THINY))
#print(paste("THINT=",THINT))

 ## unpack input parameters passed in grids
    source("unpackListedParams.R")

 ## set up grid parameters     
    source("setupGridParams.R")

 ## get lists of relative pixel locations (row,col,t) for all data and prediction pixels, and make data locatinos into standard table
    PixelLocationObj<-makeFullPIxelLocationLists()
    footprintXYT<-cbind(PixelLocationObj$xd,PixelLocationObj$yd,PixelLocationObj$td)
    yseq<-PixelLocationObj$yseq ## this is needed if we are thinning the footprint - defines which rows in each column are included
    xseq<-PixelLocationObj$xseq
    tseq<-PixelLocationObj$tseq    
      
 ## use pixel locations to obtain absolute anisotropic distance matrices
    DistMatObj<-calculateDistanceMatrices(PixelLocationObj)

 ## populate covariance matrices
    CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
    rm(DistMatObj)
    
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
          SubsetInterimMatrices.Obj<-getSubsetInterimMatrices(cDtoDmodified.list[[FP]],cDtoPmodified.list[[FP]],cPtoP,returnL=T,NuniqueFPs,FP)
          PostMeanInterim.LISTS[[FP]]<-SubsetInterimMatrices.Obj$PostMeanInterim.list
          PostVar.LISTS[[FP]]<-SubsetInterimMatrices.Obj$PostVar.list
      }
    }

 ## initialise list of length MonthDepth+1 to hold output matrices used in algorithm (and initialise first element as blank matrix to hold current month)
    OutMATlist<-vector("list",MonthDepth+1)
    OutMATlist<-lapply(OutMATlist,length) #this initialises each element of the list as a number (zero) rather than a NULL element - necessary later
    OutMATlist[[1]]<-matrix(0,nrow=Nrows,ncol=Ncols)
    #print('On creation:')
    #print(OutMATlist)
    
 ## bundle up objects to return in series of nested lists
    footprintObj<-list("footprintLUT"=footprintLUT,"footprintXYT"=footprintXYT,"modifiedFootprintTable"=modifiedFootprintTable,"yseq"=yseq,"xseq"=xseq,"tseq"=tseq)
    rm(footprintLUT); rm(footprintXYT); rm(modifiedFootprintTable)
    KrigMatObj<-list("cPtoP"=cPtoP,"PostMeanInterim.FULLlist"=PostMeanInterim.FULLlist,"PostVar.FULLlist"=PostVar.FULLlist,"PostMeanInterim.LISTS"=PostMeanInterim.LISTS,"PostVar.LISTS"=PostVar.LISTS)
    rm(cPtoP);rm(PostMeanInterim.FULLlist);rm(PostVar.FULLlist);rm(PostMeanInterim.LISTS);rm(PostVar.LISTS)
    preLoopObj<-list("footprintObj"=footprintObj,"KrigMatObj"=KrigMatObj,"OutMATlist"=OutMATlist)
    rm(footprintObj);rm(KrigMatObj)
    # print('********HERE***************')
    # for (i in length(PostMeanInterim.LISTS[[2]])){
    #     print(class(PostMeanInterim.LISTS[[2]][[i]]))
    # }
    
 ## provide summary of list structure for later checking
    listSummary<-returnListSummary(preLoopObj,paste("listSummary_preLoopObj_original_",startRel,"_",endRel,".txt",sep=""))   
    listSummary<-returnListSummary(OutMATlist,paste("listSummary_OutMATlist_original_",startRel,"_",endRel,".txt",sep=""))      

    print(gc(TRUE))

 ## return this list
    return(preLoopObj)
}
