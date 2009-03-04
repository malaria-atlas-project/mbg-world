source("transectVariog.R")
require("geoR")

 ## loop through Nomnths months and import grids
    grids<-array(NA,c(Ncols,Nrows,Nmonths))
    for(i in 1:Nmonths){
       grids[,,i]<-as.matrix(read.table(paste("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m",i,".txt",sep="")))
    }

 ##  TEMPORAL TRANSECT VARIOGRAMS

  # plot 'true' variogram based on direct MVRNORM of transect
    MVNvarObj<-makeMVNvariogramTEMPORAL(PLOT=F,Nreal=500)
    #x11();
    plot(MVNvarObj$dist,MVNvarObj$varmodel,type="n",ylim=c(0,5),xlim=c(0,Nmonths/24),col=2,lwd=3)
    title(main=paste("T; THINX=",THINX,";THINY=",THINY,";THINT=",THINT,"\nSUBSETCOL =",SUBSETCOL,";MonthDepth = ",MonthDepth))

  # plot transect variograms
    gammaVec<-0
    cnt<-0
    Nsamples<-500
    pixelX<-sample(1:Nrows,Nsamples,replace=T)
    pixelY<-sample(1:Nrows,Nsamples,replace=T)
    for(TRANS in 1:Nsamples){
       varioObj<-transectVariog(grids[pixelX[TRANS],pixelY[TRANS],],1,ORIENT="T",PLOT=F,LINES=F)
       gammaVec<-gammaVec+varioObj$v
       cnt<-cnt+1
    }
    gammaVec<-gammaVec/cnt 
 
  # replot model 
    lines(MVNvarObj$dist,MVNvarObj$varmodel,col=2,lwd=1,lty=2)
    lines(varioObj$u,gammaVec,col=4,type="o")    
    lines(MVNvarObj$u,MVNvarObj$gammaVec,col=3,type="o") 

  # also plot time series of each sample
    #x11();
    plot(seq(1/12,Nmonths/12,by=1/12),grids[pixelX[TRANS],pixelY[TRANS],],ylim=c(-15,15),type="n",ylab="f",xlab="time")
    for(TRANS in 1:10){
       lines(seq(1/12,Nmonths/12,by=1/12),grids[pixelX[TRANS],pixelY[TRANS],],type="l")
    }
  
 

 ## HORIZONTAL TRANSECT VARIOGRAMS

  # plot 'true' variogram based on direct MVRNORM of transect
    MVNvarObj<-makeMVNvariogramVERTICAL(PLOT=F,Nreal=5)
    plot(MVNvarObj$dist,MVNvarObj$varmodel,type="n",ylim=c(0,15),xlim=c(0,Ncols*3),col=2,lwd=3)
    title(main=paste("H; THINX=",THINX,";THINY=",THINY,";THINT=",THINT,"\nSUBSETCOL =",SUBSETCOL,";ColDepth = ",ColDepth))

  # plot transect variograms
    gammaVec<-0
    cnt<-0
    for(MONTH in seq(1,Nmonths,by=2)){
       for(TRANS in seq(1,min(c(Ncols,5)),by=1)){
           varioObj<-transectVariog(grids[,,MONTH],TRANS,ORIENT="H",PLOT=F,LINES=F)
           gammaVec<-gammaVec+varioObj$v
           cnt<-cnt+1
        }     
    }
    gammaVec<-gammaVec/cnt 

  # replot model 
    lines(MVNvarObj$dist,MVNvarObj$varmodel,col=2,lwd=1,lty=2)
    lines(varioObj$u,gammaVec,col=4,type="o")    
    lines(MVNvarObj$u,MVNvarObj$gammaVec,col=3,type="o") 

 ## VERTICAL TRANSECTS 

  # plot model
  # plot 'true' variogram based on duirect MVRNORM of transect
    MVNvarObj<-makeMVNvariogramVERTICAL(PLOT=F,Nreal=5)
    plot(MVNvarObj$dist,MVNvarObj$varmodel,type="n",ylim=c(0,15),xlim=c(0,Nrows*3),col=2,lwd=3)
    title(main=paste("V; THINX=",THINX,";THINY=",THINY,";THINT=",THINT,"\nSUBSETCOL =",SUBSETCOL,";ColDepth = ",ColDepth))

  # plot transect variograms
    gammaVec<-0
    cnt<-0
    for(TRANS in seq(1,min(c(Nrows,5)),by=1)){
        for(MONTH in seq(1,Nmonths,by=2)){
           varioObj<-transectVariog(grids[,,MONTH],TRANS,ORIENT="V",PLOT=F,LINES=F)
           gammaVec<-gammaVec+varioObj$v
           cnt<-cnt+1
        }
    }
    gammaVec<-gammaVec/cnt 
    
  # replot model 
    lines(MVNvarObj$dist,MVNvarObj$varmodel,col=2,lwd=1,lty=2)
    lines(varioObj$u,gammaVec,col=4,type="o") 
    lines(MVNvarObj$u,MVNvarObj$gammaVec,col=3,type="o") 






    
