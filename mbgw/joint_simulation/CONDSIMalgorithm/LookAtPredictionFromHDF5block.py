# example call
# run LookAtPredictionFromHDF5block None None None None "/home/pwg/Realizations/160779_AMS1_directsim_krigin_withlowrankerrorflag_inAndOut_spreadrel/realizations_mem_100000000_QRYPFPR220708_Americas_Run_1.9.2008_iterations_7_8.hdf5" True True True True

#run LookAtPredictionFromHDF5block None None None None "/home/pwg/Realizations/finalRealizationsForMethodsPaper/af/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_100_101.hdf5" True True True NULL "/home/pwg/Realizations/finalRealizationsForMethodsPaper/af/AMrel1.asc"


# import libraries
import numpy as np
import tables as tb
import pylab as pl
import pymc as pm
import sys
from rpy import *
from mbgw import correction_factors

# set some parameters
N_facs = int(1e5)
a_lo = 2
a_hi  = 10

# deal with system arguments
startRel = None
endRel = None
startMonth = None
endMonth = None
BACKTRANSFORM = False
AGECORRECT = False
ADDNUGGET = False

if sys.argv[1]!='None': startRel = int(sys.argv[1])
if sys.argv[2]!='None': endRel = int(sys.argv[2])
if sys.argv[3]!='None': startMonth = int(sys.argv[3])
if sys.argv[4]!='None': endMonth = int(sys.argv[4])
filename = sys.argv[5]
if sys.argv[6] == 'True': BACKTRANSFORM=True
if sys.argv[7] == 'True': AGECORRECT=True
if sys.argv[8] == 'True': ADDNUGGET=True
PLOTTING = sys.argv[9]
EXPORTASCII = sys.argv[10]
lim5kmbnry_path = sys.argv[11]

# get realization block
hf = tb.openFile(filename)
hr = hf.root

# deal with dimensions (including any slicing arguments to realizations and months)
n_realizations_RAW = hr.realizations.shape[0]
n_months_RAW = hr.realizations.shape[3]
n_rows = hr.realizations.shape[1]
n_cols = hr.realizations.shape[2]

if(startRel is None): startRel = 0
if(endRel is None): endRel = n_realizations_RAW+1
if(startMonth is None): startMonth = 0
if(endMonth is None): endMonth = n_months_RAW+1

n_realizations = (endRel-startRel)-1
n_months = (endMonth-startMonth)-1


# initialise empty block to house realisatoins of back-transformed time-aggregated PR predictions
annualmean_block = np.empty(n_rows*n_cols*n_realizations).reshape(n_rows,n_cols,n_realizations)

# optionally get nugget variance and age-correction factors    
if ADDNUGGET is True: V = hr.PyMCsamples.col('V')[:]    
if AGECORRECT is True: facs = correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)    

# loop through each realisation
for ii in xrange(0,n_realizations):

    #print ii

    # Pull out relevent section of hdf5 f block
    tot_slice = (slice(ii,ii+1,None),slice(None,None,None),slice(None,None,None),slice(startMonth,endMonth,None))  
    chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
    subsetmonth=0 
    for mm in xrange(n_months):
        chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
        subsetmonth=subsetmonth+1
    chunk = chunk.squeeze()

    holdshape = chunk.shape
    chunk = chunk.ravel()
    
    # optionally, add nugget, inverse logit, and age correct
    if ADDNUGGET is True: chunk = chunk + np.random.normal(loc=0, scale=np.sqrt(V[ii]), size=np.prod(chunk.shape))
    if BACKTRANSFORM is True: chunk = pm.invlogit(chunk)
    if AGECORRECT is True: chunk *= facs[np.random.randint(N_facs, size=np.prod(chunk.shape))]

    chunk = chunk.reshape(holdshape).squeeze()
   
    # aggregate through time
    chunkTMEAN = np.atleast_2d(np.mean(chunk,-1))
        
    # add this realisation to output block
    annualmean_block[:,:,ii]=chunkTMEAN

# get posterior mean and std of predicted maps
annualmean_mean = np.atleast_2d(np.mean(annualmean_block,-1))
annualmean_std = np.atleast_2d(np.std(annualmean_block,-1))

print 'surface mean of annual mean is '+str(np.mean(annualmean_mean))

# optionally plot
if PLOTTING == "PYLAB":

    pl.figure(1)
    pl.clf()
    pl.imshow(annualmean_mean,interpolation='nearest')
    pl.colorbar()
    pl.title('PR (posterior mean)')
    pl.axis('image')

    pl.figure(2)
    pl.clf()
    pl.imshow(annualmean_std,interpolation='nearest')
    pl.colorbar()
    pl.title('PR (posterior std)')
    pl.axis('image')
    
    pl.figure(3)
    pl.clf()
    pl.hist(annualmean_block,bins=25)
    pl.title('time-aggregated PR per-pixel (all realizations)\nmean = '+str(np.mean(annualmean_block)))

    aa = np.zeros(12).reshape(2,6)
    aa[0,1] = 1
    pl.figure(4)
    pl.clf()
    pl.imshow(aa,interpolation='nearest')
    pl.colorbar()
    pl.title('dot at top row, 2nd column ??')
    pl.axis('image')
   
    pl.show()

if PLOTTING == "RPLOT":

    # import R function
    r.source('~/mbg-world/mbgw-scripts/extract_Rlib.R')
    plotMapPY=r['plotMap']
    annualmean_mean[0,0]<-0
    annualmean_mean[0,1]<-1    
    
    template = np.zeros(1850*1850).reshape(1850,1850)
    template[0:annualmean_mean.shape[0]:1,0:annualmean_mean.shape[1]:1]=annualmean_mean
    
    plotMapPY(1-template,flipVertical="FALSE",NODATA=1,titsuf="")

if EXPORTASCII is not "NULL":

    from map_utils import getAsciiheaderFromTemplateHDF5
    from map_utils import exportAscii

    mask = tb.openFile(lim5kmbnry_path)
    hdrDict = getAsciiheaderFromTemplateHDF5(lim5kmbnry_path)

    # export as ascii
    exportAscii(annualmean_mean,EXPORTASCII,hdrDict,mask = mask.root.data[:,:])
    


    

    #r.require('fields')
    #r.image_plot(annualmean_mean)

    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()
    

#r.segments(0.5,0.05,0.5,0.95)
#r.segments(0.4,0.05,0.4,0.95)
#r.segments(np.repeat(0.4,11),np.array([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]),np.repeat(0.5,11),np.array([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]))


    
    

