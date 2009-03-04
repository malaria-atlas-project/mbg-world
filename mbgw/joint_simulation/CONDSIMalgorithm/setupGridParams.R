
 # define Ncols directly
   Ncols<-NCOLS

 # define rows and subsets
   if(class(SUBSETCOL)=="NULL") SUBSETCOL<-NROWS  # if SUBSETCOL is NULL, reset to become equal to number of rows (i.e. no subsetting of column)
   NSUBSETStemp<-NROWS/SUBSETCOL 
   NSUBSETS<-ceiling(NSUBSETStemp)
   Nrows<-NSUBSETS*SUBSETCOL

 # did we need to expand grid to make an integer number of subsets in each column?
   # print(is.null(Nrows))
   # print(Nrows)
   #    print(is.null(NROWS))
   if((floor(NSUBSETStemp)-NSUBSETStemp)!=0 & VERBOSE>=1) print(paste("Expanded grid by",Nrows-NROWS,"rows to fit subsets")) 

   if(VERBOSE>=1){
     print(paste("Imput grid: NROWS =",NROWS,"; NCOLS =",NCOLS,"; YLLCORNER =",YLLCORNER,"; CELLSIZE =",CELLSIZE))
     print(paste("Output grid: rows = ",Nrows,"; cols =",Ncols,sep=" "))
     print(paste("There will be",NSUBSETS,"subsets of",SUBSETCOL,"pixels per column"))
     print("Finished grid set-up")
    }
    
 ## calculate top edge position
    TOPEDGELAT<-YLLCORNER + (CELLSIZE*(Nrows))
