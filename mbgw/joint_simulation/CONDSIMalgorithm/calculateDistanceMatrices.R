## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################

calculateDistanceMatrices<-function(PixelLocationObj){

if(VERBOSE>=1)print("Populating distance matrices")

 ## define vectors
    xp<-PixelLocationObj$xp
    yp<-PixelLocationObj$yp
    tp<-PixelLocationObj$tp
    xd<-PixelLocationObj$xd
    yd<-PixelLocationObj$yd
    td<-PixelLocationObj$td

 ## convert pixel positions first to lat (absolute location because DDs change with lat) and long (relative location) , then to radians
    xd<-xd * CELLSIZE * (pi/180) 
    yd<-(TOPEDGELAT - (yd * CELLSIZE) + (0.5*CELLSIZE)) * (pi/180)
#    yd<-yd * CELLSIZE * (pi/180) 
    xp<-xp * CELLSIZE * (pi/180)
    yp<-(TOPEDGELAT - (yp * CELLSIZE) + (0.5*CELLSIZE)) * (pi/180)

 ## convert month vectors into appropriate temporal unit (decimal years?)
    td<-td/TEMPORALUNIT 
    tp<-tp/TEMPORALUNIT

 ## define number of data and prediciton pixels
    ndata<-length(xd) 
    npred<-length(xp) 
     
 ## populate spatial distance matrices   

  # populate D-D spatial distance matrix DtoD.iso (input/output in radians, ouptut then converted to kms)      
    geographic.list<-.Fortran("geographic",      
                  D=as.double(rep(0,ndata^2)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)),      
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata),      
                  symm=as.logical(TRUE))      
    DtoD.iso<-geographic.list$D*6378.137 # convert from radians to kms using Earth's radius      
    if(VERBOSE>=1) print("done geographic DtoD")
  
  # populate D-D bearing matrix DtoD.a               
    euc.angle.list<-.Fortran("euc_angle",      
                  theta=as.double(rep(0,ndata^2)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)), 
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata), 
                  symm=as.logical(TRUE))      
    DtoD.a<-euc.angle.list$theta   
    if(VERBOSE>=1) print("done bearing DtoD")
        
  # calculate resulting D-D anisotropic great-circle distance matrix slagDtoD     
    dist.to.aniso.list<-.Fortran("dist_to_aniso",         
                  out=as.double(rep(0,ndata^2)),         
                  D=as.double(DtoD.iso),         
                  theta=as.double(DtoD.a),         
                  nx=as.integer(ndata),         
                  ny=as.integer(ndata),         
                  inc=as.double(inc),         
                  ecc=as.double(ecc),         
                  symm=as.logical(TRUE))         
    slagDtoD<-dist.to.aniso.list$out/Scale            
    if(VERBOSE>=1) print("done aniso DtoD")
    rm(DtoD.iso)
    rm(DtoD.a)

  # populate D-P spatial distance matrix/vector DtoP.iso      
    geographic.list<-.Fortran("geographic",      
                  D=as.double(rep(0,ndata*npred)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xp,yp)),      
                  nx=as.integer(ndata),      
                  ny=as.integer(npred),      
                  symm=as.logical(FALSE))      
    DtoP.iso<-geographic.list$D*6378.137  # convert from radians to kms using Earth's radius           
    if(VERBOSE>=1) print("done geographic DtoP")

  # populate D-P bearing matrix DtoP.a     
    euc.angle.list<-.Fortran("euc_angle",      
                  theta=as.double(rep(0,ndata*npred)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xp,yp)),      
                  nx=as.integer(ndata),      
                  ny=as.integer(npred),       
                  symm=as.logical(FALSE))      
    DtoP.a<-euc.angle.list$theta      
    if(VERBOSE>=1) print("done bearing DtoP")

  # calculate resulting D-P anisotropic great-circle distance matrix slagDtoP 
    dist.to.aniso.list<-.Fortran("dist_to_aniso",
                  out=as.double(rep(0,ndata*npred)),
                  D=as.double(DtoP.iso),
                  theta=as.double(DtoP.a),
                  nx=as.integer(ndata),      
                  ny=as.integer(npred), 
                  inc=as.double(inc),
                  ecc=as.double(ecc),
                  symm=as.logical(FALSE))
    slagDtoP<-dist.to.aniso.list$out/Scale 
    if(VERBOSE>=1) print("done aniso DtoP")
    rm(DtoP.iso)
    rm(DtoP.a)

  # populate P-P spatial distance matrix PtoP.iso      
    geographic.list<-.Fortran("geographic",      
                  D=as.double(rep(0,npred*npred)),      
                  x=as.double(cbind(xp,yp)),      
                  y=as.double(cbind(xp,yp)),      
                  nx=as.integer(npred),      
                  ny=as.integer(npred),      
                  symm=as.logical(TRUE))      
    PtoP.iso<-geographic.list$D*6378.137 # convert from radians to kms using Earth's radius      
    if(VERBOSE>=1) print("done geographic PtoP")
  
  # populate P-P bearing matrix PtoP.a               
    euc.angle.list<-.Fortran("euc_angle",      
                  theta=as.double(rep(0,npred*npred)),      
                  x=as.double(cbind(xp,yp)),      
                  y=as.double(cbind(xp,yp)),      
                  nx=as.integer(npred),      
                  ny=as.integer(npred), 
                  symm=as.logical(TRUE))      
    PtoP.a<-euc.angle.list$theta      
    if(VERBOSE>=1) print("done bearing PtoP")

  # calculate resulting P-P anisotropic great-circle distance matrix slagPtoP 
    dist.to.aniso.list<-.Fortran("dist_to_aniso",
                   out=as.double(rep(0,npred*npred)),
                   D=as.double(PtoP.iso),
                   theta=as.double(PtoP.a),
                   nx=as.integer(npred),      
                   ny=as.integer(npred), 
                   inc=as.double(inc),
                   ecc=as.double(ecc),
                   symm=as.logical(TRUE))
    slagPtoP<-dist.to.aniso.list$out/Scale
    if(VERBOSE>=1) print("done aniso DtoP")
    rm(PtoP.iso)
    rm(PtoP.a)
    
 ## calculate temporal distance matrices
    Tdd<-matrix(td,nrow=length(td),ncol=length(td))
    Tpp<-matrix(tp,nrow=length(tp),ncol=length(tp))
    Tdp<-matrix(td,nrow=length(td),ncol=length(tp))
    Tpd<-matrix(tp,nrow=length(tp),ncol=length(td))
        
    tlagDtoD<-as.vector(abs(Tdd-t(Tdd)))
    tlagDtoP<-as.vector(abs(Tdp-t(Tpd)))
    tlagPtoP<-as.vector(abs(Tpp-t(Tpp)))


    Obj<-list("slagDtoD"=slagDtoD,"slagDtoP"=slagDtoP,"slagPtoP"=slagPtoP,"tlagDtoD"=tlagDtoD,"tlagDtoP"=tlagDtoP,"tlagPtoP"=tlagPtoP)
    return(Obj)
}


 


