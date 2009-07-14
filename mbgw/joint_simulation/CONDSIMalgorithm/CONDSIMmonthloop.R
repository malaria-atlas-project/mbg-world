## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################


#CONDSIMmonthloop<-function(month,preLoopObj,OutMATlist,startRel,endRel){

CONDSIMmonthloop<-function(month,preLoopObj,OutMATlist, startRel,endRel,paramfileINDEX){

    preLoopObj <- as.list(preLoopObj)
    OutMATlist <- as.list(OutMATlist)


 ## check list consistency: preLoopObj
    ListSummaryORIGINAL<-read.table(paste("listSummary_preLoopObj_original_",startRel,"_",endRel,".txt",sep=""),header=F)
    listSummary<-returnListSummary(preLoopObj,paste("listSummary_preLoopObj_new_",startRel,"_",endRel,".txt",sep="")) 
    ListSummaryNEW<-read.table(paste("listSummary_preLoopObj_new_",startRel,"_",endRel,".txt",sep=""),header=F)
    TEST<-identical(ListSummaryORIGINAL,ListSummaryNEW) 
    if(!TEST){
      print("!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!") 
      print("preLoopObj IS NOT CONSISTENT")       
      print("!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!")       
    }

 ## check list consistency: OutMATlist
    ListSummaryORIGINAL<-read.table(paste("listSummary_OutMATlist_original_",startRel,"_",endRel,".txt",sep=""),header=F)
    listSummary<-returnListSummary(OutMATlist,paste("listSummary_OutMATlist_new_",startRel,"_",endRel,".txt",sep="")) 
    ListSummaryNEW<-read.table(paste("listSummary_OutMATlist_new_",startRel,"_",endRel,".txt",sep=""),header=F)
    TEST<-identical(ListSummaryORIGINAL,ListSummaryNEW) 
    if(!TEST){
      print("!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!") 
      print("OutMATlist IS NOT CONSISTENT")       
      print("!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!")    
    }
        


 ## define inputs from bundled list
    footprintLUT<-preLoopObj$footprintObj$footprintLUT
    footprintXYT<-preLoopObj$footprintObj$footprintXYT
    modifiedFootprintTable<-preLoopObj$footprintObj$modifiedFootprintTable
    xseq<-preLoopObj$footprintObj$xseq
    yseq<-preLoopObj$footprintObj$yseq
    tseq<-preLoopObj$footprintObj$tseq

    cPtoP<-preLoopObj$KrigMatObj$cPtoP
    PostMeanInterim.FULLlist<-preLoopObj$KrigMatObj$PostMeanInterim.FULLlist
    PostVar.FULLlist<-preLoopObj$KrigMatObj$PostVar.FULLlist
    PostMeanInterim.LISTS<-preLoopObj$KrigMatObj$PostMeanInterim.LISTS
    PostVar.LISTS<-preLoopObj$KrigMatObj$PostVar.LISTS
    # for (l in PostMeanInterim.LISTS){
    #     print(class(l))
    # }


 ## check modifiedFootprintTable is logical, if not convert
    if(!is.logical(modifiedFootprintTable)){
      temp<-modifiedFootprintTable
      temp<-matrix(T,nrow=nrow(modifiedFootprintTable),ncol=ncol(modifiedFootprintTable))
      temp[modifiedFootprintTable==0]<-F
      modifiedFootprintTable<-temp
      rm(temp) 
    }

 ## invoke parameter file
#    source("ParamFile_uncond.R")
    source(paste("ParamFile_uncond_",paramfileINDEX,".R",sep=""))
    
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
         # print(month)

         # if the standard footprint:
         if(footprintType==0){ footprint<-footprintXYT}

         # otherwise make extractions from standard footprint using modifiedFootprintTable
         if(footprintType!=0){ footprint<-footprintXYT[modifiedFootprintTable[,footprintType],]}
 
         ## extract data for this col/month given footprint type
         dataMonths<-unique(footprint[,3])

         for(i in dataMonths){
            dataCols<-unique(col+footprint[footprint[,3]==i,1])
            dataRows<-yseq
            if (month==2){
              # print('*************************************')
              # print(OutMATlist[[(-i)+1]])
              # print(dataRows)  
              # print(dataCols)
#              print(col)
#              print(footprint[footprint[,3]==i,1])
            } 
            Zdata<-c(Zdata,as.vector(OutMATlist[[(-i)+1]][dataRows,dataCols]))
         } 

         ## identify correct list of interim kriging matrices for this footprint type
         # if the standard footprint:
         if(footprintType==0){
           PostMeanInterim.list<-PostMeanInterim.FULLlist
           # print('case 0')
           # print(postMeanInterim.list)           
           PostVar.list<-PostVar.FULLlist
         }     

         # if not the standard footprint:
         if(footprintType!=0){
           # print('case positive')
           # print(footprintType)     
           # print(length(PostMeanInterim.LISTS))
           # print(nrow(PostMeanInterim.LISTS[[footprintType]]))      
           # print(ncol(PostMeanInterim.LISTS[[footprintType]]))                 
           PostMeanInterim.list<-PostMeanInterim.LISTS[[footprintType]]
           # print(nrow(PostMeanInterim.list))           
           # print(ncol(PostMeanInterim.list))
           # print(class(PostMeanInterim.list))                      
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

   if(any(is.nan(MonthGrid))) print("ERROR!!! in CONDSIMmonthloop, found nan values in MonthGrid") 
   if(any(is.na(MonthGrid))) print("ERROR!!! in CONDSIMmonthloop, found na values in MonthGrid") 
   if(any(is.infinite(MonthGrid))) print("ERROR!!! in CONDSIMmonthloop, found Inf values in MonthGrid")       



    
 ## this current month now becomes month t-1, other MonthDepth months all shuffle along
  # loop through matrices in memory and shuffle them all back in time by one (if we are storing any previous months)
    if(MonthDepth>0){ 
      for(i in MonthDepth:1){
#print(paste('on i=',i,' MonthDepth=',MonthDepth,' length(OutMATlist)=',length(OutMATlist),sep="")) 
         OutMATlist[[i+1]]<-OutMATlist[[i]]
      }
    }

 ## initialise a new blank matrix for the next current month
    OutMATlist[[1]]<-matrix(0,nrow=Nrows,ncol=Ncols)
    # print('On return:')
    # print(OutMATlist)
 
 ## return grid for this month and OutMATlist object to pass on to subsequent months
    monthObject<-list("MonthGrid"=MonthGrid,"OutMATlist"=OutMATlist)
    listSummary<-returnListSummary(OutMATlist,paste("listSummary_OutMATlist_original_",startRel,"_",endRel,".txt",sep="")) 
    # image(MonthGrid)
    
    return(monthObject)
}    

