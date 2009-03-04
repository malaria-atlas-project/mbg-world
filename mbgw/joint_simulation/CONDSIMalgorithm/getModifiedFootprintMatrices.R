## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################

getModifiedFootprintMatrices<-function(footprintLUTObj,CovMatObj){

 ## define footprint LUTs
    modifiedFootprintTable<-footprintLUTObj$modifiedFootprintTable  ## table of unique modified footprints (as logical vectors of pixels in standard footprint to keep)
    NuniqueFPs<-footprintLUTObj$NuniqueFPs  ## how many unique footprints are we dealing with?
    
 ## define full covariance matrices
    cPtoP<-CovMatObj$cPtoP
    cDtoD<-CovMatObj$cDtoD
    cDtoP<-CovMatObj$cDtoP

 ## initialise list of modified interim matrices (with blank first element- special case of first col month 1)
    cDtoDmodified.list<-list(c()) 
    cDtoPmodified.list<-list(c())

 ## for each unique footprint type (but ignoring the very first - first col of month 1 is special case of no data)...
    for(i in 2:NuniqueFPs){

       ## subset matrices using footprint logical 'keep' vector
          keep<-footprintLUTObj$modifiedFootprintTable[,i]
          cDtoDmodified.list[[i]]<-cDtoD[keep,keep]
          cDtoPmodified.list[[i]]<-cDtoP[keep,]
    }
 
 ## return lists of interim output
    Obj<-list("cDtoDmodified.list"=cDtoDmodified.list,"cDtoPmodified.list"=cDtoPmodified.list)
}


