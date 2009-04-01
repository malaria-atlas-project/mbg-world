# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################

# import python libraries
import numpy as np
import copy as cp
import pymc as pm
import mbgw
import tables as tb
import time
from math import sqrt

from mbgw.joint_simulation import *
import st_cov_fun
from mbgw import correction_factors
from map_utils import getAsciiheaderFromTemplateHDF5
from map_utils import exportAscii

# import parameters from param file
from extract_params import *

# check filepaths stated in parameter file
from map_utils import checkAndBuildPaths
checkAndBuildPaths(filename,VERBOSE=True,BUILD=False)
checkAndBuildPaths(PPoutputPath,VERBOSE=True,BUILD=True)
checkAndBuildPaths(lim5kmbnry_path,VERBOSE=True,BUILD=True)

#############################################################################################################################################

#slices=[slice(None,None,None), slice(None,None,None), slice(0,12,None)]
#a_lo=2
#a_hi=10
#n_per=1
#startRel=0
#endRel=1

def extractPerPixelSummaries (slices,a_lo,a_hi,n_per,startRel,endRel):
    
    # define basic parameters
    slices = tuple(slices)     
    hf = tb.openFile(filename)    
    hr = hf.root    
    #n_realizations = len(hr.realizations)    
    n_realizations = (endRel - startRel)
    n_rows=len(hr.lat_axis)
    n_cols=len(hr.lon_axis)
    N_facs = int(1e5)
    
    N = float(n_realizations * n_per)

    # Get nugget variance and age-correction factors    
    V = hr.PyMCsamples.col('V')[:]    
    facs = mbgw.correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)    

    # import limits mask to neaten up final arrays, and check dimensions match realisation block
    mask = tb.openFile(lim5kmbnry_path)
    # perform check that the number of rows and columns is the same in both 1km grids    
    if len(mask.root.lat) != n_rows:    
        print 'WARNING!! 1km row numbers do not correspond: mask has '+str(len(mask.root.lat))+' and realisation block has '+str(n_rows)    
    if len(mask.root.long) != n_cols:    
        print 'WARNING!! col numbers do not correspond: mask has '+str(len(mask.root.long))+' and realisation block has '+str(n_cols)    

    # define a blank array of zeroes of same size as a single monthly map - that will be duplicated for various uses later
    zeroMap = np.zeros(n_rows*n_cols).reshape(n_rows,n_cols)

    # initialise zero matrices that will house running totals
    sumPR = cp.deepcopy(zeroMap)
    sumPR2 = cp.deepcopy(zeroMap)
    
    # initialise dictionary to house probability of class membership running arrays for each sceme/class
    Nschemes=len(breaksDict)    
    schemeNames=breaksDict.keys()    
    PCMdict=cp.deepcopy(breaksDict)
    
    ## ..loop through each classification scheme 
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = PCMdict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames) 

        # define additional sub-dictionary to add to PCMdict to house arrays for PCM per class per scheme
        PCM = {}

        # .. for each class within each scheme..
        for cc in xrange(0,Nclasses):
            thisbreakname = breaknames[cc]
            
            # define an empty array for this scheme-class to house PCM
            blankarray = {thisbreakname: cp.deepcopy(zeroMap) }

            # add this blank array to interim PAR dictionary
            PCM.update(blankarray)

        # add this sub-dictionary to PCMdict for this scheme
        PCM = {'PCM':PCM}
        PCMdict[scheme].update(PCM)

    # loop through each realisation
    for ii in xrange(0,n_realizations): #1:500 realisations n_realizations   
    
        # define which realisation this relates to in global set from MCMC
        MCMCrel = startRel+ii 

        print "realisation "+str(MCMCrel)+" of set "+str(startRel)+" to "+str(endRel)+" ("+str(n_realizations)+" realisations)"    

        #xxx3a = r.Sys_time() 
        # Pull out parasite rate chunk (i.e. import n months of block)    
        tot_slice = (slice(MCMCrel,MCMCrel+1,None),) + slices    
        #f_chunk = hr.realizations[tot_slice][::-1,:,:,:].T   # hr.realizations=[rel,row,col,month]   #f_chunk = [month,col,row,rel]
        #f_chunk = f_chunk[:,:,:,0]                       #f_chunk = [month,col,row]

        # because pyTables is not working properly, manually loop through each month we want to extract and make a ST 3d matrix
        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_cols,n_rows,n_months)
        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,mm] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
        f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0]   



        ########TEMP###########
        #set missing vlaues in f block to 0
        #print(sum(np.isnan(f_chunk)))
        f_chunk[np.isnan(f_chunk)]=0
        #print(sum(np.isnan(f_chunk)))
        ####################################
        
        # run check that there are no missing values in this f chunk
        if sum(sum(sum(np.isnan(f_chunk))))>0:
            print "WARNING!! found "+str(sum(np.isnan(f_chunk)))+" NaN's in realisation "+str(MCMCrel)+" EXITING!!!"
            return(-9999)

        # loop through n_per draws of the nugget..
        cnt=0
        for kk in xrange(0,n_per):

            cnt=cnt+1
            if cnt==1:
                print '    on nug rel '+str(kk)+' of '+str(n_per)
                cnt=0 
            
            # add nugget component, apply inverse logit, apply age-correction factor
            #xxx8a = r.Sys_time()
            chunk = f_chunk + np.random.normal(loc=0, scale=np.sqrt(V[MCMCrel]), size=f_chunk.shape)
            chunk = pm.invlogit(chunk.ravel())
            chunk *= facs[np.random.randint(N_facs, size=np.prod(chunk.shape))]
            chunk = chunk.reshape(f_chunk.shape).squeeze()
            #xxx8 = xxx8 + (r.Sys_time() - xxx8a)

            # aggregate through time to obtain spatial-only array for this nugget-realisation
            chunkTMEAN = np.atleast_2d(np.mean(chunk,0))
            #print('sum of chunkTMEAN '+str(sum(sum(chunkTMEAN))))


            #hdrDict = getAsciiheaderFromTemplateHDF5(lim5kmbnry_path)
            #exportAscii(chunkTMEAN,PPoutputPath+"chunkTmean.asc",hdrDict)
            #return()
            
            # increment runing total matrices 
            sumPR = sumPR + chunkTMEAN
            sumPR2 = sumPR2 + np.square(chunkTMEAN)
            
            # increment running class membership total arrays for each scheme/class
            ## ..loop through each classification scheme 
            for ss in xrange(0,Nschemes): 
                scheme=schemeNames[ss]   
                breaknames = PCMdict[scheme]['BREAKNAMES']
                Nclasses=len(breaknames) 
                breaks = PCMdict[scheme]['BREAKS']
               
                # .. for each class within each scheme..
                for cc in xrange(0,Nclasses):
                    thisbreakname = breaknames[cc]
        
                    # define an ID matrix to identify those pixels in this map in this class
                    classID = cp.deepcopy(zeroMap)
                    classID[(chunkTMEAN>=breaks[cc]) & (chunkTMEAN < breaks[cc+1])]=1 

                    # update class membership total array
                    PCMdict[scheme]['PCM'][thisbreakname] = PCMdict[scheme]['PCM'][thisbreakname] + (classID/N)

    # derive array of posterior mean
    meanPR = sumPR/N
    meanPR2 = sumPR2/N
    
    # derive array of posterior std dev
    varPR = meanPR2 - np.square(meanPR)
    stdevPR = np.sqrt(varPR)
    
    # export arrays as asciis

    ## first need to define header parameters, so copy these from an appropriate existing ascii in hdf5 format
    hdrDict = getAsciiheaderFromTemplateHDF5(lim5kmbnry_path)

    ## export meanPR as ascii
    exportAscii(meanPR,PPoutputPath+"meanPR.asc",hdrDict,mask = mask.root.data[:,:])

    ## export stdevPR as ascii
    exportAscii(stdevPR,PPoutputPath+"stdevPR.asc",hdrDict,mask = mask.root.data[:,:])
    
    ## for each classification scheme, define an array showing PCM to most likely class (PCMMLC) and what that most likely class is (MLC)
    for ss in xrange(0,Nschemes):             
        scheme=schemeNames[ss]                
        breaknames = PCMdict[scheme]['BREAKNAMES']             
        Nclasses=len(breaknames)
        
        # initialise arrays of PCMMLC and MLC            
        PCMMLC = cp.deepcopy(zeroMap)
        MLC = cp.deepcopy(zeroMap)
        
        # .. for each class within each scheme..            
        for cc in xrange(0,Nclasses):            
            thisbreakname = breaknames[cc]

            # update MLC if this class has higher PCM than previous highest
            MLCid = PCMdict[scheme]['PCM'][thisbreakname]>PCMMLC
            MLC[MLCid] = cc+1

            # keep running maximum PCM through the classes
            PCMMLC = np.maximum(PCMMLC,PCMdict[scheme]['PCM'][thisbreakname])

            # whilst at this loop location, export PCM for this scheme/class as ascii
            exportAscii(PCMdict[scheme]['PCM'][thisbreakname],PPoutputPath+'PCM_'+scheme+'_'+thisbreakname+'.asc',hdrDict,mask = mask.root.data[:,:])

        # export MLC and PCMMLC for this scheme as asciis
        exportAscii(PCMMLC,PPoutputPath+'PCMMLC_'+scheme+'.asc',hdrDict,mask = mask.root.data[:,:])
        exportAscii(MLC,PPoutputPath+'MLC_'+scheme+'.asc',hdrDict,mask = mask.root.data[:,:])
        
    return()  

#############################################################################################################################################

extractPerPixelSummaries ([slice(None,None,None), slice(None,None,None), slice(0,12,None)],2,10,10,0,100)