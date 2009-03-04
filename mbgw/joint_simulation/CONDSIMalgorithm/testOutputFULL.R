source("transectVariog.R")
require("geoR")
##############################################################################################################################
calculateVariogram<-function(Xlocs,Ylocs,data,DIR=c(),TOL=c(),PLOT=F,LINES=F,MAXDIST){


 ## calculate position of data X coordinates
     xd<-Xlocs*CELLSIZE * (pi/180)
     yd<-rep(0,length(xd))
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
      Dx<-dist.to.aniso.list$out  
      Dx<-matrix(Dx,nrow=ndata,ncol=ndata)

 ## calculate position of data Y coordinates
     yd<-(YLLCORNER + (CELLSIZE*Nrows)) - (Ylocs*CELLSIZE) * (pi/180)  #yd<-(YLLCORNER + (CELLSIZE*Nrows)) - (Ylocs*CELLSIZE) * (pi/180)
     xd<-rep(0,length(yd))
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
      Dy<-dist.to.aniso.list$out  
      Dy<-matrix(Dy,nrow=ndata,ncol=ndata)



      obj<-cbind(Dx[,1],Dy[,1],data)
      obj<-obj[!(duplicated(as.matrix(cbind(Dx[,1],Dy[,1])))),]
      geoObj<-as.geodata(obj)
      if(class(DIR)=="NULL")varioObj<-variog(geoObj,max.dist=MAXDIST,uvec=100)#breaks=seq(0,1000,length=30))
      if(class(DIR)!="NULL")varioObj<-variog(geoObj,max.dist=MAXDIST,uvec=100,direction=DIR,tolerance=TOL,unit.angle="degrees")

      
      if(PLOT) plot(varioObj$u,varioObj$v,type="o",ylim=c(0,max(varioObj$v)))
      if(LINES) lines(varioObj$u,varioObj$v,type="o") 
      return(varioObj)
} 
#####################################################################################################################################################

checkVgram<-function(monthgrid,name,MAXLAG){
    

 ## maximum lag in pixels to calculate for    
   
    gammaY<-rep(NA,MAXLAG)
    gammaX<-rep(NA,MAXLAG)

  ## define trimmed subsets that can be moved around full grid  
     subgridX<-monthgrid[,1:(ncol(monthgrid)-MAXLAG)]
     subgridY<-monthgrid[1:(nrow(monthgrid)-MAXLAG),]

  ## caclulate difference for all pixels at each lag position, and calculate gamma    
     for(i in 1:MAXLAG){
print(paste("on lag",i,"of",MAXLAG))

        startRow<-i+1
        endRow<-(nrow(monthgrid)-MAXLAG)+i 
        startCol<-i+1
        endCol<-(ncol(monthgrid)-MAXLAG)+i 

        diffGridX<-subgridX - monthgrid[,startCol:endCol] 
        diffGridY<-subgridY - monthgrid[startRow:endRow,]         

        gammaX[i]<-0.5*mean(diffGridX^2)
        gammaY[i]<-0.5*mean(diffGridY^2)        
    
    } 
    
  ## calculate corresponding theoretical models
     xp<-rep(1,MAXLAG)
     yp<-1:MAXLAG
     tp<-rep(1,MAXLAG)
     PixelLocationObj<-list("xp"=xp,"yp"=yp,"tp"=tp,"xd"=0,"yd"=0,"td"=0)
     DistMatObj<-calculateDistanceMatrices(PixelLocationObj)
     CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
     cPtoP<-CovMatObj$cPtoP
     varmodelY<-amp^2 - cPtoP[1,]
     distY<-DistMatObj$slagPtoP[1:sqrt(length(DistMatObj$slagPtoP))]*Scale

     yp<-rep(1,MAXLAG)
     xp<-1:MAXLAG
     tp<-rep(1,MAXLAG)
     PixelLocationObj<-list("xp"=xp,"yp"=yp,"tp"=tp,"xd"=0,"yd"=0,"td"=0)
     DistMatObj<-calculateDistanceMatrices(PixelLocationObj)
     CovMatObj<-PopulateCovarianceMatrices(DistMatObj)
     cPtoP<-CovMatObj$cPtoP
     varmodelX<-amp^2 - cPtoP[1,]
     distX<-DistMatObj$slagPtoP[1:sqrt(length(DistMatObj$slagPtoP))]*Scale    
   
    x11();plot(distX,gammaX,ylim=c(0,max(c(gammaX,gammaY))),type="l")
    lines(distX,varmodelX,lty=2)
    lines(distY,gammaY,type="l",col=3)
    lines(distY,varmodelY,lty=2,col=3) 
    title(main=name)   
    
}    

#####################################################################################################################################################


   ## import simulated months
    month1<-read.table("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m1.txt",header=F)
    month2<-read.table("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m2.txt",header=F)
    month3<-read.table("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m3.txt",header=F)
    month4<-read.table("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m4.txt",header=F)

 ## plot simulated months
    x11(width=15,height=15)
    par(mfrow=c(2,2))
    image(t(month1[nrow(month1):1,][1:100,1:100]),main=paste("month1"))
    image(t(month2[nrow(month1):1,][1:100,1:100]),main=paste("month2"))
    image(t(month3[nrow(month1):1,][1:100,1:100]),main=paste("month3"))
    image(t(month4[nrow(month1):1,][1:100,1:100]),main=paste("month4"))
    
   
month1smooth2<-month1
Nsmooth<-2
for(row in (Nsmooth+1):nrow(month1smooth2)){
print(row)
    month1smooth2[row,]<-colMeans(month1[row:(row-Nsmooth),])
}
month1smooth22<-month1smooth2  
for(col in (Nsmooth+1):ncol(month1smooth2)){
print(col)
    month1smooth22[,col]<-rowMeans(month1smooth2[,col:(col-Nsmooth)])
}  
    


   par(mfrow=c(1,3))
   image(t(month1[nrow(month1):1,]),main=paste("month1"))
   image(t(month1smooth22[nrow(month1):1,]),main=paste("month1smooth22"))

checkVgram(month1smooth22,"month1smooth22",200)
checkVgram(month1,"month1",200)
checkVgram(month2,"month2",200)
checkVgram(month3,"month3",200)
checkVgram(month4,"month4",200)
  
  