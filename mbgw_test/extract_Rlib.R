########################################################################################################## 
plotMap<-function(inputmatrix,titsuf="x",NODATA=c()){
  if(class(inputmatrix)!="matrix"){
    print(paste("WARNING!! input to PlotMap is not a matrix (",class(inputmatrix),") - will not be plotted",sep=""))
    return()
  }
  if(class(NODATA)!="NULL") inputmatrix[inputmatrix==NODATA]<-NA
  x11();image(t(inputmatrix[nrow(inputmatrix):1,]),main=paste("month",titsuf))
}
########################################################################################################## 
plotmonth<-function(grid,dummy=c(),dummy2=c()){
 
 if(class(grid)=="array"){
   Nmonths=dim(grid)[1]
   print(paste("Nmonths",Nmonths))
   
   for(t in 1:Nmonths){
      gridmonth=grid[t,,]
      plotMap(gridmonth,t)
   }
 }
 if(class(grid)=="matrix") plotMap(grid)
}
##########################################################################################################  
getTimeMean<-function(inputarray,TimeDim){
  
 #print(class(inputarray))
 #print(dim(inputarray))
 
 # create time-averaged grid accross all time slices of passed cube
 meanGrid<-colMeans(inputarray,dims=TimeDim)
 
 #naindex<-is.na(meanGrid)
 #if(sum(naindex!=0)){
 #  print(paste("WARNING!! Found",sum(naindex),"NAs in time-averaged matrix"))
 #  #meanGrid[naindex]<- -9999
 #}

 return(meanGrid)     
} 
########################################################################################################## 
flipVertical<-function(inputmatrix){

  if(class(inputmatrix)!="matrix"){
    print(paste("WARNING!! input to flipVertical is not a matrix (",class(inputmatrix),") - returned unchanged",sep=""))
    return()
  }

  outputmatrix<-inputmatrix[nrow(inputmatrix):1,]  

  return(outputmatrix)
}
##########################################################################################################
expandGridRes<-function(inputarray,expFac){

  #print(inputarray)
  #print(class(inputarray))
  #print(dim(inputarray))

  if(class(inputarray)=="numeric"){
    inputarray=array(inputarray)
  }

  dims<-dim(inputarray)
  Nd<-length(dims)
  ncolIN=dims[Nd]

  if(length(dims)>2){
    print(paste("WARNING!!! input to expandGridRes has >2 dimensions (",length(dims),") EXITING",sep=""))
    return()
  }

  if(Nd==1){ 
    colIndex<-as.vector(t(matrix(1:ncolIN,nrow=ncolIN,ncol=expFac)))
    outputarray=inputarray[colIndex]
  }
  
  if(Nd==2){
    nrowIN=dims[1]
    colIndex<-as.vector(t(matrix(1:ncolIN,nrow=ncolIN,ncol=expFac)))
    rowIndex<-as.vector(t(matrix(1:nrowIN,nrow=nrowIN,ncol=expFac)))
    outputarray=inputarray[rowIndex,colIndex]
  }

  return(outputarray)
}
##########################################################################################################

