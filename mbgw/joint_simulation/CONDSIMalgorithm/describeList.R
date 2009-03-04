## Author: Dr Peter Gething, University of Oxford
## Date: 6th February 2009
## License: Creative Commons BY-NC-SA
##################################################

#############################################################################################
getListStructure<-function(LIST,level){
  Nelements<-length(LIST)
  elementnames<-names(LIST)
  if(class(elementnames)!="NULL")namesOrder<-order(elementnames)
  if(class(elementnames)=="NULL")namesOrder<-1:Nelements  

  for(i in namesOrder){
     if(class(LIST[[i]])=="list"){
       if(class(elementnames)!="NULL")rawTable<<-rbind(rawTable,c(level,length(LIST[[i]]),"list",elementnames[[i]]))
       if(class(elementnames)=="NULL")rawTable<<-rbind(rawTable,c(level,length(LIST[[i]]),"list","NONE"))
       getListStructure(LIST[[i]],level+1)
     }
     if(class(LIST[[i]])!="list"){
        CLASS<-"unknown" 
        if(class(LIST[[i]])=="NULL") CLASS<-"NULL"        
        if(is.matrix(LIST[[i]])) CLASS<-"matrix"
        if(is.vector(LIST[[i]])) CLASS<-"vector"
        if(is.vector(LIST[[i]]) & length(LIST[[i]])==1) CLASS<-"scalar"
       if(class(elementnames)!="NULL")rawTable<<-rbind(rawTable,c(level,0,CLASS,elementnames[[i]]))
       if(class(elementnames)=="NULL")rawTable<<-rbind(rawTable,c(level,0,CLASS,"NONE"))
     }
  }
} 
#############################################################################################
returnListSummary<-function(INPUTLIST,OUTPUTPATH=c()){

  rawTable<<-NULL 
  getListStructure(INPUTLIST,1)
  rawTable<-rbind(c(0,length(INPUTLIST),"list","NONE"),rawTable)

  formattedTable<-as.data.frame(matrix(NA,nrow=nrow(rawTable),ncol=ncol(rawTable)+max(as.numeric(rawTable[,1])-1)))
  for(i in 1:nrow(rawTable)){
     newRow<-c(rep("",rawTable[i,1]),rawTable[i,-1])
     widthRow<-length(newRow)
     formattedTable[i,1:widthRow]<-newRow

  }

  formattedTable[is.na(formattedTable)]<-""

  if(class(OUTPUTPATH)!="NULL") write.table(formattedTable,OUTPUTPATH,row.names=F,col.names=F,quote=F)
  return(formattedTable)
} 
#############################################################################################
  