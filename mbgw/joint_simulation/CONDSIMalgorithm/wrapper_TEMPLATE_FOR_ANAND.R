
## 13 SCALAR VALUES YOU NEED TO DEFINE...
Scale<-510.4378
amp<-3.289207
inc<-1.179354
ecc<-0.3091333
t.lim.corr<-0.6484932
scale.t<-4.183328
sin.frac<-0.03048946

YLLCORNER<--33.08455
CELLSIZE<-0.04166665
NROWS<-1481 #100
NCOLS<-1752 #100
 
Nmonths<-5
StartMonth<-1985.0 - 2009


## ..THEN BUNDLE INTO THREE LISTS
covParamObj<-list("Scale"=Scale,"amp"=amp,"inc"=inc,"ecc"=ecc,"t.lim.corr"=t.lim.corr,"scale.t"=scale.t,"sin.frac"=sin.frac)
gridParamObj<-list("YLLCORNER"=YLLCORNER,"CELLSIZE"=CELLSIZE,"NROWS"=NROWS,"NCOLS"=NCOLS)
monthParamObj<-list("Nmonths"=Nmonths,"StartMonth"=StartMonth)


####################################################################################################################################################################


## THEN ALL OF THE BELOW NEEDS TO BE REPLICATED IN PYTHON..

 ## source two main functions 
    source("CONDSIMpreloop.R")
    source("CONDSIMmonthloop.R")

 ## call pre-loop function to define necessary matrices etc
    preLoopObj<-CONDSIMpreloop(covParamObj,gridParamObj,monthParamObj)

  # initialise OutMATlist from pre-loop object
    OutMATlist<-preLoopObj$OutMATlist

 ## main loop through Nmonth months
    for(month in 1:Nmonths){

       ## call main CONDSIMmonthloop function
       monthObject<-CONDSIMmonthloop(month,preLoopObj,OutMATlist)
       OutMATlist<-monthObject$OutMATlist
       MonthGrid<-monthObject$MonthGrid
    }
    
## SO 'MonthGrid' IS YOUR UNCONDITIONAL REALISATION OF THE GRID FOR MONTH month...TO DO WITH AS YOU WISH ;-)    

