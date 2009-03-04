source("transectVariog.R")
require("geoR")

 ## loop through Nomnths months and import grids
    grids<-array(NA,c(Ncols,Nrows,Nmonths))
    for(i in 1:Nmonths){
       grids[,,i]<-as.matrix(read.table(paste("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m",i,".txt",sep="")))
    }

 ##  TRANSECT VARIOGRAMS

  # plot 'true' variogram based on direct MVRNORM of transect
    MVNvarObj<-makeMVNvariogramTEMPORAL(PLOT=F,Nreal=500)
    #x11();
    plot(MVNvarObj$dist,MVNvarObj$varmodel,type="n",ylim=c(0,5),xlim=c(0,Nmonths/24),col=2,lwd=3)
    title(main=paste("T; THIN =",THIN,";SUBSETCOL =",SUBSETCOL,";ColDepth =",ColDepth,";MonthDepth =",MonthDepth))

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
  
 