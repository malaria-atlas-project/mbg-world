from rpy import *
import numpy as np
import copy as cp
#import pymc as pm
import tables as tb
#import mbgw


def getGridCovarianceInY(gridIN,meanIN,cellWidth,VERBOSE=False,NaNvalue="nan"):

    # ensure all NA values are set to nan
    if NaNvalue!="nan":
        gridIN[np.isnan(gridIN)]=NaNvalue
    
    # define rows to evalulate at
    n_rows = shape(gridIN)[0]
    n_cols = shape(gridIN)[1]

    YpixList = np.array(xrange(0,n_rows))    
    littleJumps = YpixList[0:int(n_rows*0.25):2]
    bigJumps  = YpixList[int(n_rows*0.25):int(n_rows*0.8):20]
    jumps=np.concatenate((littleJumps,bigJumps))
    #jumps=YpixList[0:int(n_rows):2]

    # initialise vector to house covariance at each lag
    E_cov = np.empty(len(jumps)-1)
    
    # define interim grid
    tempGrid =  gridIN-meanIN

    # loop through pixel jumps and compare two offset grids at that displacement
    for d in xrange(len(jumps)-1):
   
        if(VERBOSE==True): print 'on jump '+str(d)+' of '+str(len(jumps))
        bumper=np.zeros((jumps[d],n_cols))
        bumper[bumper==0]=nan
    
        grid1 =  np.vstack((tempGrid,bumper))    
        grid2 =  np.vstack((bumper,tempGrid))
        tempProd = grid1*grid2
        E_cov[d] = np.mean(tempProd[np.isnan(tempProd)==False])
        
    # calculate distance in radians at the displacement locations
    RadDist = jumps[:-1:]*cellWidth
    
    # return dictionary
    retDict = {"PixDist":jumps[:-1:],"RadDist":RadDist,"E_cov":E_cov}
    return (retDict)
