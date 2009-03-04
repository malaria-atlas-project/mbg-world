## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################

getPixelSubset<-function(PixelLocationObj,col,month){

 ## check if we need to trim any pixel locations from the data footprint, checking respectively:
    # not full set of columns to left (current or previous months) OR
    # not full set of columns to right (previous months only) OR  
    # not full set of previous months  
    if(!((col<=ColDepth) | ((MonthDepth>0 & month>1) & ((Ncols-col) < ColDepth)) | (month<=MonthDepth))) return()
    
 ## define vectors
    xd<-PixelLocationObj$xd
    yd<-PixelLocationObj$yd
    td<-PixelLocationObj$td

 ## initialise  logical vector: true = keep
    keep<-rep(T,length(xd))

 ## carry out tests on pixel location lists, marking those pixels that need to be trimmed from the footprint
 
  # footprint pixels too far to left 
    keep[xd< -(col-1)]<-F
    
  # footprint pixels too far to right (only applies on previous month layers because curent month takes only from left) 
    keep[xd>(Ncols-col)]<-F    

  # footprint pixels from a month that precedes the earliest available
    keep[td< (-(month-1))]<-F 
    
    return(keep)
}
 
