
 ## define time parameters
    TEMPORALUNIT<-12 

 ## define search parameters
    ColDepth<-1   # how many columns of previously simulated data will be used to predict subsequent column
    MonthDepth<-1 # how many months of previously simulated data will be used to predict subsequent months 
    THINX<-1
    THINY<-50
    THINT<-1
    

 ## are we going to simulate subsets of each column? non null value specifies length of each subset
    #SUBSETCOL<-50     
    SUBSETCOL<-c()

 ## do we want to print row/col progress, warnings etc
    VERBOSE<-1 #0 = no messages; 1= some messages; 2=loads of messages

 ## define parameter nu (~ matern degree of differentiability) [NB only valid if doing ST]
    nu<-0.5

 ## define whether we are definnig the covaiance function with a sinusoidal seasonal component
   SEASONAL<-TRUE 
 
