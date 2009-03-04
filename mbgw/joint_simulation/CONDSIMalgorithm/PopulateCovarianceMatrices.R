## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## Liscence: Creative Commons BY-NC-SA
##################################################

PopulateCovarianceMatrices<-function(DistMatObj){

if(VERBOSE>=1)print("Populating covariance matrices")

  ## define distance matrices locally
     slagDtoD<-DistMatObj$slagDtoD
     slagDtoP<-DistMatObj$slagDtoP
     slagPtoP<-DistMatObj$slagPtoP
     tlagDtoD<-DistMatObj$tlagDtoD
     tlagDtoP<-DistMatObj$tlagDtoP
     tlagPtoP<-DistMatObj$tlagPtoP

 ## define number of dat            
    ndata<-sqrt(length(slagDtoD))
    npred<-sqrt(length(slagPtoP))
    
#  # data to data
#    cDtoD<-exp(-slagDtoD)*(amp^2)
#    cDtoD<-matrix(cDtoD,nrow=ndata,ncol=ndata)      
#
#    cDtoP<-exp(-slagDtoP)*(amp^2)
#    cDtoP<-matrix(cDtoP,nrow=ndata,ncol=npred)
#
#    cPtoP<-exp(-slagPtoP)*(amp^2)
#    cPtoP<-matrix(cPtoP,nrow=npred,ncol=npred)

  ## first define temporal 'variogram' vectors or matrices               
     if(!SEASONAL){              
       GtDtoD<- 1/((t.lim.corr/nu) + (((1-t.lim.corr)/nu)*exp(-(tlagDtoD/scale.t))) )              
       GtDtoP<- 1/((t.lim.corr/nu) + (((1-t.lim.corr)/nu)*exp(-(tlagDtoP/scale.t))) )              
       GtPtoP<- 1/((t.lim.corr/nu) + (((1-t.lim.corr)/nu)*exp(-(tlagPtoP/scale.t))) )              
     }              
                   
     if(SEASONAL){              
       GtDtoD<- 1/(((1-t.lim.corr)/nu) * (sin.frac*cos(2*pi*tlagDtoD) + (1-sin.frac) * exp(-tlagDtoD /scale.t)) + (t.lim.corr/nu))              
       GtDtoP<- 1/(((1-t.lim.corr)/nu) * (sin.frac*cos(2*pi*tlagDtoP) + (1-sin.frac) * exp(-tlagDtoP /scale.t)) + (t.lim.corr/nu))              
       GtPtoP<- 1/(((1-t.lim.corr)/nu) * (sin.frac*cos(2*pi*tlagPtoP) + (1-sin.frac) * exp(-tlagPtoP /scale.t)) + (t.lim.corr/nu))              
     }                             

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
     # print(ndata)
     # print(npred)
     # print(1*cDtoD)
     # print((amp[1]^2)*matrix(c(1,0,0,1),ncol=2,nrow=2))
     # print(class(amp[1]))
     # print(mean(cDtoD))              
     cDtoD<-cDtoD*(amp^2)              
     if(VERBOSE>=1) print("done DtoD")
                  
   # data to prediction              
     stein.list<-.Fortran("stein_spatiotemporal",              
                     C=as.double(slagDtoP),              
                     Gt=as.double(GtDtoP),              
                     origin_val=as.double(nu),              
                     nx=as.integer(ndata),              
                     ny=as.integer(npred),              
                     symm=as.logical(FALSE))              
     cDtoP<-matrix(stein.list$C,nrow=ndata,ncol=npred)              
     cDtoP<-cDtoP*(amp^2)              
     if(VERBOSE>=1) print("done DtoP")
               
   # prediction to prediction              
     stein.list<-.Fortran("stein_spatiotemporal",              
                     C=as.double(slagPtoP),              
                     Gt=as.double(GtPtoP),              
                     origin_val=as.double(nu),              
                     nx=as.integer(npred),              
                     ny=as.integer(npred),              
                     symm=as.logical(TRUE))              
     cPtoP<-matrix(stein.list$C,nrow=npred,ncol=npred)              
     cPtoP<-cPtoP*(amp^2)              
     if(VERBOSE>=1) print("done PtoP")                    

   
 ## return covariance matrices            
    return(list("cDtoD"=cDtoD,"cDtoP"=cDtoP,"cPtoP"=cPtoP)) 

}
