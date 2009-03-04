## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################

getfootprintLUTs<-function(PixelLocationObj){

 ## initialise table where each col is a unique modified footprint (logical vector detialing pixels to discard from standard footprint)
    modifiedFootprintTable<-NULL
              
 ## initialise LUT matrix of Ncols by Nmonths: each cell will contain number of required footprint for this column/month
    footprintLUT<-matrix(0,ncol=Ncols,nrow=Nmonths)
  
 ## initialise counter for unique modified footprints
    uniqueFP<-0   
    
 ## loop through all columns and months to be predicted..
      
    for(month in 1:Nmonths){  #month<-1
       for(col in 1:Ncols){ #col<-2
             
           ## get data pixel subset: if we are near an edge we may not be using the full data footprint 
              pixelSubset<-getPixelSubset(PixelLocationObj,col,month)
                
           ## if modified footprint is needed (non-NULL return)..
              if(class(pixelSubset)!="NULL"){
 
                ## if first subset to be defined, start modifiedFootprintTable and increment unique footprint ID, and this mnth/col location assigned footprint ID 1.
                   if(class(modifiedFootprintTable)=="NULL"){
                     modifiedFootprintTable<-matrix(pixelSubset)
                     uniqueFP<-uniqueFP+1
                     footprintLUT[month,col]<-1
                   }

                ## if not the first subset:
                   else{
                     
                     # check if this footprint is unque, or if it has been defined before
                       matched<-apply(modifiedFootprintTable,2,identical,pixelSubset)

                     # if this is not a new unique footprint, then this mnth/col location is assigned this existing footprint ID
                       if(sum(matched)>0){
                         footprintLUT[month,col]<-which(matched)
                       }

 
                     # if this is a new unique footprint, then new footprint appended to modifiedFootprintTable, unique footprint ID incremented, and this mnth/col location assigned this new footprint ID
                       if(sum(matched)==0){
                         modifiedFootprintTable<-cbind(modifiedFootprintTable,pixelSubset)
                         uniqueFP<-uniqueFP+1
                         footprintLUT[month,col]<-uniqueFP
                       }                  

                   }
              } 
       } 
    }

    Obj<-list("footprintLUT"=footprintLUT,"modifiedFootprintTable"=modifiedFootprintTable,"NuniqueFPs"=uniqueFP)     
}
