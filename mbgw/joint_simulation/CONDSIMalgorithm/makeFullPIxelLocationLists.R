############################################################################# 
makeFullPIxelLocationLists<-function(){

 ## define vectors of prediction pixel locations 
    xp<-rep(0,Nrows)
    yp<-1:Nrows
    tp<-rep(0,Nrows)

 ## define sequence of pixel row and month positions
    yseq<-1:Nrows

  # if we are thinning the footprint, thin the yseq vector accordingly
    if(class(THINY)!="NULL"){
      yseq<-yseq[seq(1,length(yseq),by=THINY)]
    }
    
    if(MonthDepth==0) tseq<-NULL
    if(MonthDepth!=0){
      tseq<--1:-MonthDepth
      if(class(THINT)!="NULL")tseq<-tseq[seq(1,length(tseq),by=THINT)]
    }  
      

  # define lists of positions for all pixels: contemporary month
    xseq<--ColDepth:-1 # only taking ColDepth columns to left of prediction column
    if(class(THINX)=="numeric") xseq<-rev(rev(xseq)[seq(1,length(xseq),by=THINX)])   ## implement thinning at specified level - this line also ensures that, whatever thinning is imlemented, the nearst left column is always included.
    if(THINX=="EXPO") xseq<-rev(-(1:floor(sqrt(ColDepth)))^2) ## implement thinning by spacing pixels exponentially along width of ColDepth

print("xseq 1:")
print(xseq)
 
    yd<-rep(yseq,length(xseq))
    xd<-as.vector(t(matrix(xseq,nrow=length(xseq),ncol=length(yseq))))
    td<-rep(0,length(xseq)*length(yseq))

  # define lists of positions for all pixels: previous months up to MonthDepth
    if(MonthDepth!=0){
      xseq<-c(xseq,0,rev(abs(xseq)))  #taking ColDepth columns to both left and right of prediction column, plus the one directly underneath
print("xseq 2:")
print(xseq)
      yd<-c(yd,rep(rep(yseq,length(xseq)),length(tseq)))
      xd<-c(xd,rep(as.vector(t(matrix(xseq,nrow=length(xseq),ncol=length(yseq)))),length(tseq)))
      td<-c(td,as.vector(t(matrix(tseq,nrow=length(tseq),ncol=length(xseq)*length(yseq))))) 
    }

    if(VERBOSE>=1){
      print(paste("full footprint has",length(xd),"data point"))
      print(paste("predicting at",length(xp),"locations"))      
    }
    
  ## return pixel location list
     Obj<-list("xp"=xp,"yp"=yp,"tp"=tp,"xd"=xd,"yd"=yd,"td"=td,"yseq"=yseq,"xseq"=xseq,"tseq"=tseq)
     return(Obj)
}
#############################################################################
  
