from rpy import *
import numpy as np
import copy as cp
#import pymc as pm
import tables as tb
#import mbgw


def getGridCovarianceInT(cubeIN,meanIN,THINNING=20,VERBOSE=False,NaNvalue="nan"):

    # ensure all NA values are set to nan
    if NaNvalue!="nan":
        cubeIN[np.isnan(cubeIN)]=NaNvalue

    # define nshapes
    n_rows = shape(cubeIN)[0]
    n_cols = shape(cubeIN)[1]
    n_months= shape(cubeIN)[2]

    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()   


    xindex = np.arange(0,n_rows,THINNING)
    yindex = np.arange(0,n_cols,THINNING)
        
    ## loop through each month..and extract a thinned sample of pixel values, converted into a 1-darray, 
    ## which is then appended to matching vectors from every other month

    for i in arange(n_months):

        # extract a thin set of pixels from this month and convert to 1d array for convenience
        monthExtraction = cubeIN[0:n_rows:THINNING,0:n_cols:THINNING,i]
        monthExtraction=monthExtraction.ravel()
        
        # if we are dealing with a conditoined block, subtract from mean constant for this month
        #print 'len(meanIN): '+str(len(meanIN))
        if len(shape(meanIN))>=1: monthExtraction=monthExtraction-meanIN[i]       
        
        # bind to running matrix 
        if i ==0: monthMat=monthExtraction
        if i>0: monthMat=np.vstack((monthMat,monthExtraction))

    # initialise lag sizes
    jumps=np.arange(n_months)
    #jumps=YpixList[0:int(n_rows):2]

    # initialise vector to house covariance at each lag
    E_cov = np.empty(len(jumps)-1)

    # loop through pixel jumps and compare two offset grids at that displacement
    for d in xrange(len(jumps)-1):
   
        if(VERBOSE==True): print 'on jump '+str(d)+' of '+str(len(jumps))
        bumper=np.zeros((jumps[d],monthMat.shape[1]))
        bumper[bumper==0]=nan
    
        grid1 =  np.vstack((monthMat,bumper))    
        grid2 =  np.vstack((bumper,monthMat))
        tempProd = grid1*grid2
        E_cov[d] = np.mean(tempProd[np.isnan(tempProd)==False])

    # calculate distance in years at the displacement locations
    yearDist = jumps[:-1:]*(1./12)

    # return dictionary
    retDict = {"PixDist":jumps[:-1:],"yearDist":yearDist,"E_cov":E_cov}
    
    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()
    
    return (retDict)
