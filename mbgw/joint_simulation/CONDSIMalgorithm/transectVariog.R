###############################################################################################################################
transectVariog<-function(grid,TRANS,ORIENT,PLOT,LINES,COL=1){
#grid<-grids[,,MONTH]
#TRANS<-22
#ORIENT<-"H"

   if(ORIENT=="H"){
     xd<-seq(1,ncol(grid))*CELLSIZE * (pi/180)
     yd<-rep((YLLCORNER + (CELLSIZE*Nrows)) - (TRANS*CELLSIZE),ncol(grid)) * (pi/180)
     data<-grid[TRANS,]
     ndata<-length(data) 
     geographic.list<-.Fortran("geographic",      
                  D=as.double(rep(0,ndata*ndata)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)),      
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata),      
                  symm=as.logical(TRUE))      
      DtoD.iso<-geographic.list$D*6378.137 # convert from radians to kms using Earth's radius      
              
      euc.angle.list<-.Fortran("euc_angle",      
                  theta=as.double(rep(0,ndata*ndata)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)), 
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata), 
                  symm=as.logical(TRUE))      
      DtoD.a<-euc.angle.list$theta   
    
      dist.to.aniso.list<-.Fortran("dist_to_aniso",         
                  out=as.double(rep(0,ndata*ndata)),         
                  D=as.double(DtoD.iso),         
                  theta=as.double(DtoD.a),         
                  nx=as.integer(ndata),         
                  ny=as.integer(ndata),         
                  inc=as.double(inc),         
                  ecc=as.double(ecc),         
                  symm=as.logical(TRUE))         
      D<-dist.to.aniso.list$out  
      D<-matrix(D,nrow=ndata,ncol=ndata)

      obj<-cbind(D[,1],rep(0,ndata),data)
      geoObj<-as.geodata(obj)
      varioObj<-variog(geoObj,max.dist=3000,uvec=30)#breaks=seq(0,1000,length=30))
      
      if(PLOT) plot(varioObj$u,varioObj$v,type="o")
      if(LINES) lines(varioObj$u,varioObj$v,type="o",col=COL)
      return(varioObj) 
    }
 

   if(ORIENT=="V"){
     yd<-((seq(nrow(grid),1)*CELLSIZE)+YLLCORNER) * (pi/180)
     xd<-rep(TRANS,nrow(grid)) * (pi/180)
     data<-grid[,TRANS]
     ndata<-length(data) 
     geographic.list<-.Fortran("geographic",      
                  D=as.double(rep(0,ndata*ndata)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)),      
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata),      
                  symm=as.logical(TRUE))      
      DtoD.iso<-geographic.list$D*6378.137 # convert from radians to kms using Earth's radius      
              
      euc.angle.list<-.Fortran("euc_angle",      
                  theta=as.double(rep(0,ndata*ndata)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)), 
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata), 
                  symm=as.logical(TRUE))      
      DtoD.a<-euc.angle.list$theta   
    
      dist.to.aniso.list<-.Fortran("dist_to_aniso",         
                  out=as.double(rep(0,ndata*ndata)),         
                  D=as.double(DtoD.iso),         
                  theta=as.double(DtoD.a),         
                  nx=as.integer(ndata),         
                  ny=as.integer(ndata),         
                  inc=as.double(inc),         
                  ecc=as.double(ecc),         
                  symm=as.logical(TRUE))         
      D<-dist.to.aniso.list$out  
      D<-matrix(D,nrow=ndata,ncol=ndata)
      
      obj<-cbind(D[,1],rep(0,ndata),data)
      geoObj<-as.geodata(obj) 
      varioObj<-variog(geoObj,max.dist=3000,uvec=30)
      
      if(PLOT) plot(varioObj$u,varioObj$v,type="o")
      if(LINES) lines(varioObj$u,varioObj$v,type="o",col=COL)
      return(varioObj) 
    }
    
   if(ORIENT=="T"){
      D<-(1:length(grid))/12
      data<-grid
      obj<-cbind(D,rep(0,length(grid)),as.vector(data))
      geoObj<-as.geodata(obj)
      varioObj<-variog(geoObj,max.dist=Nmonths/12,breaks=seq(0,(Nmonths/12),by=1/12),uvec=30)
      
      if(PLOT) plot(varioObj$u,varioObj$v,type="o")
      if(LINES) lines(varioObj$u,varioObj$v,type="o",col=COL)
      return(varioObj) 
    }    
} 
###############################################################################################################################
makeMVNvariogramVERTICAL<-function(PLOT=F,Nreal=50){

   xp<-rep(1,Nrows)
   yp<-1:Nrows
   tp<-rep(1,Nrows)
   xd<-0
   yd<-0
   td<-0
   PixelLocationObj<-list("xp"=xp,"yp"=yp,"tp"=tp,"xd"=xd,"yd"=yd,"td"=td)
   
 
 ## use pixel locations to obtain absolute anisotropic distance matrices 
    DistMatObj<-calculateDistanceMatrices(PixelLocationObj)

 ## populate covariance matrices
    CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
    
 ## define PtoP covariance matrix (this is constant - prediction set is always a single column - regardless of where we are predicting)   
    cPtoP<-CovMatObj$cPtoP
    
 ## calculate mean vector for this prediction column
    meanP<-rep(0,nrow(cPtoP))

    grid<-matrix(NA,nrow=Nrows,ncol=Ncols) 
    
  # plot model 
    varmodel<-amp^2 - cPtoP[1,]
    dist<-DistMatObj$slagPtoP[1:sqrt(length(DistMatObj$slagPtoP))]*Scale
    
#    x11();plot(1:1000,covmodel,type="n",ylim=c(0,15),col=2,lwd=3)
    if(PLOT){
      x11();plot(dist,varmodel,type="n",ylim=c(0,15),xlim=c(0,1000),col=2,lwd=3)
      title(main="Vertical transect: MVN")
    }
      
 ## create realisations of this transect and a variogram for each 
    gammaVec<-0
    for(i in 1:Nreal){   
       grid[,1]<-MVRNORM(1,MU=meanP,COV=cPtoP)
       varioObj<-transectVariog(grid,1,ORIENT="V",FALSE,LINES=PLOT)
       gammaVec<-gammaVec+varioObj$v
    }
    gammaVec<-gammaVec/Nreal     

  # replot model  
    if(PLOT){
      lines(varioObj$u,gammaVec,col=3,type="o")
      lines(dist,varmodel,col=2,lwd=1,lty=2)
    }
    
    return(list("dist"=dist,"varmodel"=varmodel,"u"=varioObj$u,"gammaVec"=gammaVec))

}
###############################################################################################################################
makeMVNvariogramHORIZONTAL<-function(PLOT=F,Nreal=50){
   xp<-1:Ncols
   yp<-rep(1,Ncols)
   tp<-rep(1,Ncols)
   xd<-0
   yd<-0
   td<-0
   PixelLocationObj<-list("xp"=xp,"yp"=yp,"tp"=tp,"xd"=xd,"yd"=yd,"td"=td)
   
 ## use pixel locations to obtain absolute anisotropic distance matrices
    DistMatObj<-calculateDistanceMatrices(PixelLocationObj)

 ## populate covariance matrices
    CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
    
 ## define PtoP covariance matrix (this is constant - prediction set is always a single column - regardless of where we are predicting)   
    cPtoP<-CovMatObj$cPtoP
    
 ## calculate mean vector for this prediction column
    meanP<-rep(0,nrow(cPtoP))

    grid<-matrix(NA,nrow=Nrows,ncol=Ncols) 
    
  # plot model    
    varmodel<-amp^2 - cPtoP[1,]
    dist<-DistMatObj$slagPtoP[1:sqrt(length(DistMatObj$slagPtoP))]*Scale
    
#    x11();plot(1:1000,covmodel,type="n",ylim=c(0,15),col=2,lwd=3)
    if(PLOT){
      x11();plot(dist,varmodel,type="n",ylim=c(0,15),xlim=c(0,1000),col=2,lwd=3)
      title(main="Horizontal transect: MVN") 
    }   
        
 ## create realisations of this transect and a variogram for each 
    gammaVec<-0
    for(i in 1:Nreal){   
       grid[1,]<-MVRNORM(1,MU=meanP,COV=cPtoP)
       varioObj<-transectVariog(grid,1,ORIENT="H",FALSE,LINES=PLOT)
       gammaVec<-gammaVec+varioObj$v
    }
    gammaVec<-gammaVec/Nreal     

  # replot model  
    if(PLOT){
      lines(varioObj$u,gammaVec,col=3,type="o")
      lines(dist,varmodel,col=2,lwd=1,lty=2)
    }
        
    return(list("dist"=dist,"varmodel"=varmodel,"u"=varioObj$u,"gammaVec"=gammaVec))
}  
###############################################################################################################################
makeMVNvariogramTEMPORAL<-function(PLOT=F,Nreal=50){
   xp<-rep(1,Nmonths)
   yp<-rep(1,Nmonths)
   tp<-1:Nmonths
   xd<-0
   yd<-0
   td<-0
   PixelLocationObj<-list("xp"=xp,"yp"=yp,"tp"=tp,"xd"=xd,"yd"=yd,"td"=td)
   
 ## use pixel locations to obtain absolute anisotropic distance matrices
    DistMatObj<-calculateDistanceMatrices(PixelLocationObj)

 ## populate covariance matrices
    CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
    
 ## define PtoP covariance matrix (this is constant - prediction set is always a single column - regardless of where we are predicting)   
    cPtoP<-CovMatObj$cPtoP
    
 ## calculate mean vector for this prediction column
    meanP<-rep(0,Nmonths)

    grid<-rep(NA,Nmonths)
    
  # plot model     
    varmodel<-amp^2 - cPtoP[1,]
    dist<-DistMatObj$tlagPtoP[1:sqrt(length(DistMatObj$tlagPtoP))]
    
#    x11();plot(1:1000,covmodel,type="n",ylim=c(0,15),col=2,lwd=3)
    if(PLOT){
      x11();plot(dist,varmodel,type="n",ylim=c(0,15),xlim=c(0,Nmonths/12),col=2,lwd=3)
      title(main="Temporal transect: MVN") 
    }   
        
 ## create realisations of this transect and a variogram for each 
    gammaVec<-0
    for(i in 1:Nreal){   
       grid<-MVRNORM(1,MU=meanP,COV=cPtoP)
       varioObj<-transectVariog(grid,1,ORIENT="T",FALSE,LINES=PLOT)
       gammaVec<-gammaVec+varioObj$v
    }
    gammaVec<-gammaVec/Nreal      

  # replot model  
    if(PLOT){
      lines(varioObj$u,gammaVec,col=3,type="o")
      lines(dist,varmodel,col=2,lwd=1,lty=2)
    }
        
    return(list("dist"=dist,"varmodel"=varmodel,"u"=varioObj$u,"gammaVec"=gammaVec))
} 
###############################################################################################################################
