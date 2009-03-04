
##################################TEMP
Scale.DUMMY<-510.4378
amp.DUMMY<-3.289207
inc.DUMMY<-1.179354
ecc.DUMMY<-0.3091333
t.lim.corr.DUMMY<-0.6484932
scale.t.DUMMY<-4.183328
sin.frac.DUMMY<-0.03048946

YLLCORNER.DUMMY<--33.08455
CELLSIZE.DUMMY<-0.04166665
NROWS.DUMMY<-100# 1481 #
NCOLS.DUMMY<-100#1752 #100
 
Nmonths.DUMMY<-5
StartMonth.DUMMY<-1985.0 - 2009

covParamObj.DUMMY<-list("Scale"=Scale.DUMMY,"amp"=amp.DUMMY,"inc"=inc.DUMMY,"ecc"=ecc.DUMMY,"t.lim.corr"=t.lim.corr.DUMMY,"scale.t"=scale.t.DUMMY,"sin.frac"=sin.frac.DUMMY)
gridParamObj.DUMMY<-list("YLLCORNER"=YLLCORNER.DUMMY,"CELLSIZE"=CELLSIZE.DUMMY,"NROWS"=NROWS.DUMMY,"NCOLS"=NCOLS.DUMMY)
monthParamObj.DUMMY<-list("Nmonths"=Nmonths.DUMMY,"StartMonth"=StartMonth.DUMMY)
##################################TEMP

 ## source two main functions 
    source("CONDSIMpreloop.R")
    source("CONDSIMmonthloop.R")

 ## call pre-loop function to define necessary matrices etc
a<-Sys.time()
    preLoopObj<-CONDSIMpreloop(covParamObj.DUMMY,gridParamObj.DUMMY,monthParamObj.DUMMY)
A<-Sys.time()-a

  # initialise OutMATlist from pre-loop object
    OutMATlist<-preLoopObj$OutMATlist

 ## main loop through Nmonth months
b<-Sys.time()
    for(month in 1:Nmonths){
print(paste("month",month))

       ## call main CONDSIMmonthloop function
       monthObject<-CONDSIMmonthloop(month,preLoopObj,OutMATlist)
       OutMATlist<-monthObject$OutMATlist
       MonthGrid<-monthObject$MonthGrid

       ## export this current month and optionally display
       exportPath<-paste("/home/pwg/CONDSIM/CONDSIMoutput/uncongrids/uncon_k1_m",month,".txt",sep="")
       write.table(MonthGrid,exportPath,row.names=F,col.names=F)
    }
B<-Sys.time()-b

 ## print how long these processes took
    print(paste("Pre loop:",A))    
    print(paste("Loop over",Nmonths,"months:",B))    


 ## generate test variograms of output
#    source("testOutputFULL.R")
    
    


#}}}}    
