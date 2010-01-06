# example call
# run examineRealization("/home/pwg/Realizations/150709/realizations_mem_100000000_QRYPFPR220708_Americas_Run_1.9.2008_iterations_0_1.hdf5",0,11,17,None,None,False,"FALSE",True,True)

# import python libraries
from rpy import *
import numpy as np
import copy as cp
import pymc as pm
import tables as tb
import mbgw
from pr_incidence import *
import time
import st_cov_fun
from mbgw.joint_simulation import *
from math import sqrt
from mbgw import correction_factors
from map_utils import getAsciiheaderFromTemplateHDF5
from map_utils import exportAscii



from getGridCovarianceInY import * 
from getGridCovarianceInT import *
# import R function
r.source('~/mbg-world/mbgw-scripts/extract_Rlib.R')
plotMapPY=r['plotMap']
flipVerticalPY=r['flipVertical']
r.source('temptestcov.R')
temptestcovPY=r['testRcov']



def examineRealization (filename,Rel,Month,paramfileINDEX,TemporalStartMonth=None,TemporalEndMonth=None,conditioned=False,flipVertical="FALSE",SPACE=True,TIME=True):

    # deal with system arguments
    #filename = sys.argv[1]
    #Rel = int(sys.argv[2]) 
    #Month = int(sys.argv[3]) 
    #conditioned = sys.argv[4]
    #flipVertical = sys.argv[5]
    #paramfileINDEX = int(sys.argv[6])
    #TemporalStartMonth = int(sys.argv[7])
    #TemporalEndMonth = int(sys.argv[8])
    #SPACE = sys.argv[9]
    #TIME = sys.argv[10]

    ## if filename is a string, assume its a path and import the hdf5 file (otherwise, assumption is we are pasing the 'hr' root of an hdf5 realisation file)
    if type(filename) is str:
        hf = tb.openFile(filename)    
        hr = hf.root
        
    if type(filename) is not str:    
        hr = filename
        
    # define path to R param file
    mbgw_root = __root__ = mbgw.__path__[0]
    r_paramfile_path= mbgw_root+'/joint_simulation/CONDSIMalgorithm/ParamFile_uncond_'+str(paramfileINDEX)+'.R'

    # initialise plot window
    nplots = 0
    if SPACE is True: nplots=nplots+5
    if TIME is True: nplots=nplots+1
    r.X11(width=3.3*nplots,height=4)
    r.par(mfrow=(1,nplots))

    ###CHECK SPATIAL COVARIANCE AND BASIC FEATURE OF A SINGLE MONTH
    if SPACE is True:

        # define basic parameters
        slices=[slice(None,None,None), slice(None,None,None), slice(Month,Month+1,None)]

        slices = tuple(slices)     
        n_realizations = 1
        n_rows=len(hr.lat_axis)
        n_cols=len(hr.lon_axis)
        N_facs = int(1e5)

        # Pull out parasite rate chunk (i.e. import n months of block)    
        slices = tuple(slices)  
        tot_slice = (slice(Rel,Rel+1,None),) + slices    

        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
        subsetmonth=0 

        #print tot_slice
        #print f_chunk[:,:,:,subsetmonth]

        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
            subsetmonth=subsetmonth+1
        #f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0] 
        f_chunk = f_chunk.squeeze()
        f_chunk[f_chunk==-9999]=nan

        # plot this grid 
        plotMapPY(f_chunk.squeeze(),flipVertical=flipVertical)
        r.title(main="logit")

        inv_f_chunk=pm.invlogit(f_chunk.squeeze().T)
        inv_f_chunk=inv_f_chunk.reshape(shape(f_chunk))
        plotMapPY(inv_f_chunk,flipVertical=flipVertical)
        r.title(main="inv logit")


        #from IPython.Debugger import Pdb
        #Pdb(color_scheme='Linux').set_trace()   


        # compare global variance to parameter draw
        observedVar = round(np.var(f_chunk[np.isnan(f_chunk)==False]),10)
        theoreticalVar = ((hr.PyMCsamples.col('amp')[Rel])**2)
        varString = 'observedVar = :'+str(observedVar)+';  amp^2 =: '+str(theoreticalVar)
        print varString

        # plot histogram
        junk=r.hist(f_chunk[np.isnan(f_chunk)==False],main=varString,xlab="",ylab="")
        junk=r.hist(pm.invlogit(f_chunk[np.isnan(f_chunk)==False]),xlab="",ylab="", main="")

        # calculate and plot empirical covariance function in N-S direction
        gridIN = cp.deepcopy(f_chunk).squeeze()

        if conditioned is False: meanIN=0
        if conditioned is True: meanIN = hr.PyMCsamples.col("m_const")[Rel] + (hr.PyMCsamples.col("t_coef")[Rel]*hr.t_axis[Month])
        cellWidth=5/6378.137
        covDict = getGridCovarianceInY(gridIN,meanIN,cellWidth)    

        # obtain theoretical covariance function from input MCMC paramater values: pymc method
        C = hr.group0.C[Rel]
        xplot = covDict['RadDist']
        yplot1 = C([[0,0,0]], np.vstack((np.zeros(len(xplot)),xplot,np.zeros(len(xplot)))).T)
        yplot1 = np.asarray(yplot1).squeeze()

        # obtain theoretical covariance function from input MCMC paramater values: R method
        Scale=hr.PyMCsamples.col("scale")[Rel]
        amp=hr.PyMCsamples.col("amp")[Rel]
        inc=hr.PyMCsamples.col("inc")[Rel]
        ecc=hr.PyMCsamples.col("ecc")[Rel]
        t_lim_corr=hr.PyMCsamples.col("t_lim_corr")[Rel]
        scale_t=hr.PyMCsamples.col("scale_t")[Rel]
        sin_frac=hr.PyMCsamples.col("sin_frac")[Rel]

        CfromR=temptestcovPY(xplot,np.zeros(len(xplot)),np.zeros(len(xplot)),Scale,amp,inc,ecc,t_lim_corr,scale_t,sin_frac,r_paramfile_path)
        yplot = CfromR[0,:]

        # plot

        ymax = max(np.max(covDict['E_cov']),np.max(xplot),np.max(yplot))
        ymin = min(np.min(covDict['E_cov']),np.min(xplot),np.min(yplot))

        r.plot(covDict['RadDist'],covDict['E_cov'],xlab="radians",ylab="C",main=str(paramfileINDEX),ylim=(ymin,ymax))    
        r.lines(xplot,yplot1,col=2)
        r.lines(xplot,yplot,col=3)

    ###CHECK TEMPORAL COVARIANCE

    if TIME is True:

        # if start and months are None, or if they are non-valid, rest to maximum temporal extents
        if ((TemporalEndMonth is None) | (TemporalEndMonth>=hr.realizations.shape[3])): TemporalEndMonth=hr.realizations.shape[3]
        if ((TemporalStartMonth is None) | (TemporalStartMonth>=(hr.realizations.shape[3]-1))): TemporalStartMonth=0



        # define basic parameters
        slices=[slice(None,None,None), slice(None,None,None), slice(TemporalStartMonth,TemporalEndMonth,None)]

        slices = tuple(slices)     
        n_realizations = 1
        n_rows=len(hr.lat_axis)
        n_cols=len(hr.lon_axis)
        N_facs = int(1e5)


        # Pull out parasite rate chunk (i.e. import n months of block)    
        slices = tuple(slices)  
        tot_slice = (slice(Rel,Rel+1,None),) + slices    

        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
        subsetmonth=0 

        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
            subsetmonth=subsetmonth+1
        #f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0]
        f_chunk = f_chunk.squeeze()   
        f_chunk[f_chunk==-9999]=nan

        # calculate and plot empirical temporal covariance
        gridIN = cp.deepcopy(f_chunk).squeeze()

        if conditioned is False: meanIN=0
        if conditioned is True: meanIN = hr.PyMCsamples.col("m_const")[Rel] + (hr.PyMCsamples.col("t_coef")[Rel]*hr.t_axis[TemporalStartMonth:TemporalEndMonth+1:1])

        covDict = getGridCovarianceInT(gridIN,meanIN)    

        # obtain theoretical covariance function from input MCMC paramater values: pymc method
        C = hr.group0.C[Rel]
        xplot = covDict['yearDist']
        yplot = C([[0,0,0]], np.vstack((np.zeros(len(xplot)),np.zeros(len(xplot)),xplot)).T)
        yplot = np.asarray(yplot).squeeze()

        # obtain theoretical covariance function from input MCMC paramater values: R method
        Scale=hr.PyMCsamples.col("scale")[Rel]
        amp=hr.PyMCsamples.col("amp")[Rel]
        inc=hr.PyMCsamples.col("inc")[Rel]
        ecc=hr.PyMCsamples.col("ecc")[Rel]
        t_lim_corr=hr.PyMCsamples.col("t_lim_corr")[Rel]
        scale_t=hr.PyMCsamples.col("scale_t")[Rel]
        sin_frac=hr.PyMCsamples.col("sin_frac")[Rel]

        CfromR=temptestcovPY(np.zeros(len(xplot)),np.zeros(len(xplot)),xplot,Scale,amp,inc,ecc,t_lim_corr,scale_t,sin_frac,r_paramfile_path)
        yplot2 = CfromR[0,:]

        # plot

        ymax = max(np.max(covDict['E_cov']),np.max(yplot),np.max(yplot2))
        ymin = min(np.min(covDict['E_cov']),np.min(yplot),np.min(yplot2),0)

        r.plot(covDict['yearDist'],covDict['E_cov'],xlab="lag (years)",ylab="C",main=str(paramfileINDEX),ylim=(ymin,ymax))    
        r.lines(xplot,yplot,col=2)
        r.lines(xplot,yplot2,col=3)

def assessRealizationCovariance (filename,Rel,Month,paramfileINDEX,TemporalStartMonth=None,TemporalEndMonth=None,conditioned=False,flipVertical="FALSE",SPACE=True,TIME=True):

    # deal with system arguments
    #filename = sys.argv[1]
    #Rel = int(sys.argv[2]) 
    #Month = int(sys.argv[3]) 
    #conditioned = sys.argv[4]
    #flipVertical = sys.argv[5]
    #paramfileINDEX = int(sys.argv[6])
    #TemporalStartMonth = int(sys.argv[7])
    #TemporalEndMonth = int(sys.argv[8])
    #SPACE = sys.argv[9]
    #TIME = sys.argv[10]

    ## if filename is a string, assume its a path and import the hdf5 file (otherwise, assumption is we are pasing the 'hr' root of an hdf5 realisation file)
    if type(filename) is str:
        hf = tb.openFile(filename)    
        hr = hf.root
        
    if type(filename) is not str:    
        hr = filename
        
    # define path to R param file
    mbgw_root = __root__ = mbgw.__path__[0]
    r_paramfile_path= mbgw_root+'/joint_simulation/CONDSIMalgorithm/ParamFile_uncond_'+str(paramfileINDEX)+'.R'

    ###CHECK SPATIAL COVARIANCE AND BASIC FEATURE OF A SINGLE MONTH
    if SPACE is True:

        # define basic parameters
        slices=[slice(None,None,None), slice(None,None,None), slice(Month,Month+1,None)]

        slices = tuple(slices)     
        n_realizations = 1
        n_rows=len(hr.lat_axis)
        n_cols=len(hr.lon_axis)
        N_facs = int(1e5)

        # Pull out parasite rate chunk (i.e. import n months of block)    
        slices = tuple(slices)  
        tot_slice = (slice(Rel,Rel+1,None),) + slices    

        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
        subsetmonth=0 

        #print tot_slice
        #print f_chunk[:,:,:,subsetmonth]

        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
            subsetmonth=subsetmonth+1
        #f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0] 
        f_chunk = f_chunk.squeeze()
        f_chunk[f_chunk==-9999]=nan
        inv_f_chunk=pm.invlogit(f_chunk.squeeze().T)
        inv_f_chunk=inv_f_chunk.reshape(shape(f_chunk))

        #from IPython.Debugger import Pdb
        #Pdb(color_scheme='Linux').set_trace()   

        # calculate empirical covariance function in N-S direction
        gridIN = cp.deepcopy(f_chunk).squeeze()

        if conditioned is False: meanIN=0
        if conditioned is True: meanIN = hr.PyMCsamples.col("m_const")[Rel] + (hr.PyMCsamples.col("t_coef")[Rel]*hr.t_axis[Month])
        cellWidth=5/6378.137
        covDict = getGridCovarianceInY(gridIN,meanIN,cellWidth)    

        # obtain theoretical covariance function from input MCMC paramater values: pymc method
        C = hr.group0.C[Rel]
        xplot = covDict['RadDist']
        yplot1 = C([[0,0,0]], np.vstack((np.zeros(len(xplot)),xplot,np.zeros(len(xplot)))).T)
        yplot1 = np.asarray(yplot1).squeeze()

#        # obtain theoretical covariance function from input MCMC paramater values: R method
#        Scale=hr.PyMCsamples.col("scale")[Rel]
#        amp=hr.PyMCsamples.col("amp")[Rel]
#        inc=hr.PyMCsamples.col("inc")[Rel]
#        ecc=hr.PyMCsamples.col("ecc")[Rel]
#        t_lim_corr=hr.PyMCsamples.col("t_lim_corr")[Rel]
#        scale_t=hr.PyMCsamples.col("scale_t")[Rel]
#        sin_frac=hr.PyMCsamples.col("sin_frac")[Rel]

#        CfromR=temptestcovPY(xplot,np.zeros(len(xplot)),np.zeros(len(xplot)),Scale,amp,inc,ecc,t_lim_corr,scale_t,sin_frac,r_paramfile_path)
#        yplot = CfromR[0,:]

        # plot
        
        Slag_emp = covDict['RadDist']
        Slag_mod = xplot
        Scov_emp = covDict['E_cov']
        Scov_mod = yplot1
        

    ###CHECK TEMPORAL COVARIANCE

    if TIME is True:

        # if start and months are None, or if they are non-valid, rest to maximum temporal extents
        if ((TemporalEndMonth is None) | (TemporalEndMonth>=hr.realizations.shape[3])): TemporalEndMonth=hr.realizations.shape[3]
        if ((TemporalStartMonth is None) | (TemporalStartMonth>=(hr.realizations.shape[3]-1))): TemporalStartMonth=0

        # define basic parameters
        slices=[slice(None,None,None), slice(None,None,None), slice(TemporalStartMonth,TemporalEndMonth,None)]

        slices = tuple(slices)     
        n_realizations = 1
        n_rows=len(hr.lat_axis)
        n_cols=len(hr.lon_axis)
        N_facs = int(1e5)


        # Pull out parasite rate chunk (i.e. import n months of block)    
        slices = tuple(slices)  
        tot_slice = (slice(Rel,Rel+1,None),) + slices    

        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
        subsetmonth=0 

        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
            subsetmonth=subsetmonth+1
        #f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0]
        f_chunk = f_chunk.squeeze()   
        f_chunk[f_chunk==-9999]=nan

        # calculate and plot empirical temporal covariance
        gridIN = cp.deepcopy(f_chunk).squeeze()

        if conditioned is False: meanIN=0
        if conditioned is True: meanIN = hr.PyMCsamples.col("m_const")[Rel] + (hr.PyMCsamples.col("t_coef")[Rel]*hr.t_axis[TemporalStartMonth:TemporalEndMonth+1:1])

        covDict = getGridCovarianceInT(gridIN,meanIN)    

        # obtain theoretical covariance function from input MCMC paramater values: pymc method
        C = hr.group0.C[Rel]
        xplot = covDict['yearDist']
        yplot = C([[0,0,0]], np.vstack((np.zeros(len(xplot)),np.zeros(len(xplot)),xplot)).T)
        yplot = np.asarray(yplot).squeeze()

#        # obtain theoretical covariance function from input MCMC paramater values: R method
#        Scale=hr.PyMCsamples.col("scale")[Rel]
#        amp=hr.PyMCsamples.col("amp")[Rel]
#        inc=hr.PyMCsamples.col("inc")[Rel]
#        ecc=hr.PyMCsamples.col("ecc")[Rel]
#        t_lim_corr=hr.PyMCsamples.col("t_lim_corr")[Rel]
#        scale_t=hr.PyMCsamples.col("scale_t")[Rel]
#        sin_frac=hr.PyMCsamples.col("sin_frac")[Rel]

#        CfromR=temptestcovPY(np.zeros(len(xplot)),np.zeros(len(xplot)),xplot,Scale,amp,inc,ecc,t_lim_corr,scale_t,sin_frac,r_paramfile_path)
#        yplot2 = CfromR[0,:]

        # plot

        Tlag_emp = covDict['yearDist']
        Tlag_mod = xplot
        Tcov_emp = covDict['E_cov']
        Tcov_mod = yplot

        retDict = {'Slag_emp':Slag_emp,'Slag_mod':Slag_mod,'Scov_emp':Scov_emp,'Scov_mod':Scov_mod,'Tlag_emp':Tlag_emp,'Tlag_mod':Tlag_mod,'Tcov_emp':Tcov_emp,'Tcov_mod':Tcov_mod}
        return(retDict)

    retDict = {'Slag_emp':Slag_emp,'Slag_mod':Slag_mod,'Scov_emp':Scov_emp,'Scov_mod':Scov_mod}
    return(retDict)

def compileMultipleCovarianceComparisons(realizationsFolderPath):

    relCovDict=[]
    for fname in os.listdir(realizationsFolderPath):

        # is this file and hdf5 file?
        if fname.find('.hdf5') == -1: continue
        print 'on file: '+fname

        relCovDict.append(assessRealizationCovariance (realizationsFolderPath+fname,0,11,17,TemporalStartMonth=None,TemporalEndMonth=None,conditioned=False,SPACE=True,TIME=True))

    return relCovDict

def plotCovarianceComparisons(relCovDict,SPACE=True,TIME=True,toPDF=None):


    # initialise plot window
    nplots = 0
    if SPACE is True: nplots=nplots+2
    if TIME is True: nplots=nplots+2
    if toPDF is not None: r.pdf(toPDF,width=nplots*4,height=4)
    if toPDF is  None: r.X11(width=nplots*4,height=4)
    r.par(mfrow=(1,nplots))

    # plot ensemble of covariance functions (all found in realizationsFolderPath)
    if SPACE is True:
        r.plot((0,0.5),(-50,50),main="Empirical space-marginal covariance functions",ylab="covariance",xlab="lag (radians)",type="n")
        for i in xrange(len(relCovDict)):
            r.lines(relCovDict[i]['Slag_emp'],relCovDict[i]['Scov_emp'])

        r.plot((0,0.5),(-50,50),main="Theoretical space-marginal covariance functions",ylab="covariance",xlab="lag (radians)",type="n")
        for i in xrange(len(relCovDict)):
            r.lines(relCovDict[i]['Slag_mod'],relCovDict[i]['Scov_mod'])

    if TIME is True:
        r.plot((0,1),(-50,50),main="Empirical time-marginal covariance functions",ylab="covariance",xlab="lag (years)",type="n")
        for i in xrange(len(relCovDict)):
            r.lines(relCovDict[i]['Tlag_emp'],relCovDict[i]['Tcov_emp'])

        r.plot((0,1),(-50,50),main="Theoretical time-marginal covariance functions",ylab="covariance",xlab="lag (years)",type="n")
        for i in xrange(len(relCovDict)):
            r.lines(relCovDict[i]['Tlag_mod'],relCovDict[i]['Tcov_mod'])

    if toPDF is not None: r.dev_off()







