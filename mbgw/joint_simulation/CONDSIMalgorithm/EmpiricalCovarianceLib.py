import numpy as np
import copy as cp
import pymc as pm
import tables as tb
from rpy import *
import mbgw
import st_cov_fun
from mbgw.joint_simulation import *
from math import sqrt


###########TEMP
#x=np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
#y=(x*2)
#xy=np.vstack((x,y)).T
#t= 2008-(np.array([1,1,2,3,2,5,3,4,6,11])/12.)
#z = np.array([0.1,0.2,0.1,-0.4,-0.6,0,0.3,0.25,-0.4,-0.2])
#mu=0
#nbins=10
################

def getEmpiricalCovarianceFunction(xy,z,mu,nbins=None):

    # work out if we passed a two- or one-column set of locations
    if (len(np.shape(xy))==1):
        d1=xy
        d2=None
        
    if (len(np.shape(xy))==2):
        d1=xy[:,0]
        d2=xy[:,1]
        

    # get distance matrix

    ## if no second dimension supplied, assume we are evaluating along a transect, or through time:
    if (d2 is None):
        X = np.vstack(((d1,))*len(d1))
        dxy = X-X.T
        
    ## if two dimensions supplied, assum,e we want cross-distances (i.e. we have two spatial dimensions)
    if (d2 is not None):
        X = np.vstack(((d1,))*len(d1))
        dX = X-X.T
        Y = np.vstack(((d2,))*len(d2))
        dY = Y-Y.T
        Dxy = np.sqrt(dX*dX +dY*dY)

    # get empirical point-to-point covariance matrix
    temp = z-mu
    TEMP = np.vstack(((temp,))*len(temp))
    TEMP = TEMP*TEMP.T
    
    # convert lag and C matrices to vector based only on lower triangle(and diagonal)
    IDmat = np.ones(np.product(np.shape(TEMP))).reshape(np.shape(TEMP))
    Cvector = TEMP[np.where(np.tril(IDmat,k=0)==1)]
    Dvector = Dxy[np.where(np.tril(IDmat,k=0)==1)]

    # if we are not binning, then return vectors of all pairwise covariances and lags
    if (nbins is None): return ({'C':Cvector,'lag':Dvector})
    
    # if we are binning
    
    ## establish maximum distance, and get nbins equal bins along this distance
    mxlag = Dvector.max()
    binWidth = mxlag/(nbins)
    binMins = np.arange(nbins)*binWidth
    binMaxs = binMins+binWidth
    
    # loop through bins and get expected covariance
    Cbins = np.ones(nbins)*-9999
    Dbins = np.ones(nbins)*-9999
    for i in np.arange(nbins):
        lagID = np.where((Dvector>=binMins[i]) & (Dvector<binMaxs[i]))
        if (len(lagID[0])>0):  # i.e only populate this lag if there is any data in it
            Cbins[i]=np.mean(Cvector[lagID])
            Dbins[i]=np.mean(Dvector[lagID])
        
    # remove bins with no data
    IDrm = np.where((Cbins==-9999) | (Dbins==-9999))
    
    return ({'C':Cbins,'lag':Dbins})

def plotEmpiricalCovarianceFunction(EmpCovFuncDict,CovModelObj=None, cutoff = 0.8, title=None):

    # plot empirical covariance function
    if (cutoff is None): XLIM=(0,EmpCovFuncDict['lag'].max())
    if (cutoff is not None): XLIM=(0,EmpCovFuncDict['lag'].max()*cutoff)

    ymin = np.min((EmpCovFuncDict['C'].min(),0))    
    ymax = EmpCovFuncDict['C'].max()


    # if passed, add plot title
    if (title is not None):
        r.title(main=title)
        
    # if passed, add theoretical covariance model
    if (CovModelObj is not None):
        xplot = EmpCovFuncDict['lag']
        yplot = CovModelObj([[0,0,0]], np.vstack((np.zeros(len(xplot)),xplot,np.zeros(len(xplot)))).T)
        yplot = np.asarray(yplot).squeeze()
        r.lines(xplot,yplot,col=3)    
        if (np.max(yplot)>ymax): ymax = np.max(yplot)
        if (np.min(yplot)<ymin): ymin = np.min(yplot)

    YLIM = (ymin,ymax)    
    r.plot(EmpCovFuncDict['lag'],EmpCovFuncDict['C'],ylim=YLIM,xlim=XLIM,xlab="lag",ylab="covariance")

