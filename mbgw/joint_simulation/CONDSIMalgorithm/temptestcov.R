
testRcov<-function(xd,yd,td,Scale,amp,inc,ecc,t.lim.corr,scale.t,sin.frac,r_paramfile_path){

#Scale=0.0800292937636586
#amp=3.28920723639468
#inc=1.17935378082411
#ecc=0.309133346084652
#t.lim.corr=0.648493215546719
#scale.t=4.18332793066728
#sin.frac=0.0304894629404626

#xd=c(0.1, 0.2, 0.3, 0.4, 0.5)
#yd=c(0.1, 0.2, 0.3, 0.4, 0.5)
#td = c(0,0,0,0,0)



#print(paste("Scale=",Scale,sep=""))
#print(paste("amp=",amp,sep=""))
#print(paste("inc=",inc,sep=""))
#print(paste("ecc=",ecc,sep=""))
#print(paste("t.lim.corr=",t.lim.corr,sep=""))
#print(paste("scale.t=",scale.t,sep=""))
#print(paste("sin.frac=",sin.frac,sep=""))

#print(xd)
#print(yd)
#print(td)

#d<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
#d<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
#d<-rep(0,10)


#xd<-xd * CELLSIZE * (pi/180) 
#yd<-(TOPEDGELAT - (yd * CELLSIZE) + (0.5*CELLSIZE)) * (pi/180)
#td<-td/TEMPORALUNIT 

ndata<-length(xd)

#source("ParamFile_uncond.R")


source(r_paramfile_path)

 ## source fortran functions
    temp<-system("R CMD SHLIB geographic.f",intern=T)
    temp<-system("R CMD SHLIB euc_angle.f",intern=T)
    temp<-system("R CMD SHLIB dist_to_aniso.f",intern=T)
    temp<-system("R CMD SHLIB isotropic_cov_funs.f",intern=T)
    temp<-system("R CMD SHLIB Stein_st_covariance_serial.f",intern=T)
    dyn.load("geographic.so")
    dyn.load("euc_angle.so")
    dyn.load("dist_to_aniso.so")
    dyn.load("isotropic_cov_funs.so")
    dyn.load("Stein_st_covariance_serial.so") 

  # populate D-D spatial distance matrix DtoD.iso (input/output in radians, ouptut then converted to kms)      
    geographic.list<-.Fortran("geographic",      
                  D=as.double(rep(0,ndata^2)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)),      
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata),      
                  symm=as.logical(TRUE))      
    DtoD.iso<-geographic.list$D


  # populate D-D bearing matrix DtoD.a               
    euc.angle.list<-.Fortran("euc_angle",      
                  theta=as.double(rep(0,ndata^2)),      
                  x=as.double(cbind(xd,yd)),      
                  y=as.double(cbind(xd,yd)), 
                  nx=as.integer(ndata),      
                  ny=as.integer(ndata), 
                  symm=as.logical(TRUE))      
    DtoD.a<-euc.angle.list$theta   


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


 ## calculate temporal distance matrices
    Tdd<-matrix(td,nrow=length(td),ncol=length(td))
    tlagDtoD<-as.vector(abs(Tdd-t(Tdd)))



 ## define number of dat            
    ndata<-sqrt(length(slagDtoD))


  ## first define temporal 'variogram' vectors or matrices               
       GtDtoD<- 1/(((1-t.lim.corr)/nu) * (sin.frac*cos(2*pi*tlagDtoD) + (1-sin.frac) * exp(-tlagDtoD /scale.t)) + (t.lim.corr/nu))              

  ## now populate covariance vectors/matrices using stein ST covariance function               
      
   # data to data              
     stein.list<-.Fortran("stein_spatiotemporal",              
                     C=as.double(slagDtoD),              
                     Gt=as.double(GtDtoD),              
                     origin_val=as.double(nu),              
                     nx=as.integer(ndata),              
                     ny=as.integer(ndata),              
                     symm=as.logical(TRUE))              
     cDtoD<-matrix(stein.list$C,nrow=ndata,ncol=ndata)
     cDtoD<-cDtoD*(amp^2)
     
     return("cDtoD"=cDtoD)
}         

