# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################

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

rowsPerChunk=1

# import R function
r.source('extract_Rlib.R')
expandGridResPY=r['expandGridRes']

# import parameters from param file
from extract_params import *
 
#############################################################################################################################################
def examineSalb (salblim1km,uniqueSalb_path={},pixelN_path={},ignore={}):

    ''''
    Takes an input raster with integers defining unique spatial areas (e.g. the salblim1km grid)
    and calculates two vectors: uniqueSalb is a list of unique integer IDs present in this raster,
    and count is the corresponding number of pixels in each. These are returned as a dictionary
    and are optionally exported to txt files.
    
    Params to pass are:
    
    salblim1km      : either an hdf5 filepath or a numpy array object containing the salb raster
    uniqueSalb_path : filepath to export uniqueSalb - no export if ommited
    pixelN_path     : filepath to export count - no export if ommited
    ignore          : a list containing values of salb to ignore from output lists (e.g. -9999)
    '''

    # if the salb file is specified with a filepath, import it  
    if type(salblim1km) == str:
        salblim1km = tb.openFile(salblim1km, mode = "r")

    nrow=shape(salblim1km.root.data)[0]
    uniqueSalb=[-9999]
    count=[0]

    for i in xrange(0,nrow):
        #print(r.paste("row",i,"of",nrow))
        chunk=salblim1km.root.data[i,:]
        uniqueSalbCHUNK=unique(chunk)
        NuniqueSalbCHUNK=len(uniqueSalbCHUNK)
        countCHUNK= zeros(NuniqueSalbCHUNK)
        for j in xrange(0,NuniqueSalbCHUNK):

            countCHUNK[j]=sum(chunk==uniqueSalbCHUNK[j])

            if uniqueSalbCHUNK[j] in uniqueSalb:
                index=where(uniqueSalb==uniqueSalbCHUNK[j])[0]
                count[index]=count[index]+countCHUNK[j]
            else:
                uniqueSalb=append(uniqueSalb,uniqueSalbCHUNK[j]) 
                count=append(count,countCHUNK[j])

    # optionally remove stated values (e.g 0 or -9999) from list of unique salbs
    if len(ignore)>0:
        for i in xrange(0,len(ignore)):
            count = count[uniqueSalb!=ignore[i]] 
            uniqueSalb = uniqueSalb[uniqueSalb!=ignore[i]] 

    # optionally export to file
    if len(uniqueSalb_path)>0:
        uniqueSalb.tofile(uniqueSalb_path,sep=",")
    if len(pixelN_path)>0:    
        count.tofile(pixelN_path,sep=",")
        
    # return unique list and corresponding pixel count as a dictionary
    returnDict={'uniqueSalb':uniqueSalb,'count':count}
    return returnDict
#############################################################################################################################################
def extractSummaries_country(slices,a_lo,a_hi,n_per,FileStartRel,FileEndRel,startRel=None,endRel=None):

    '''
    Takes an hdf block of one or more realisation of f, and calculates mean PR, total burden, and population at risk (PAR) in each
    unique spatial unit (e.g. country) specified in the 1km res salblim1km_path file. Also requires 1km population surface specified by 
    grump1km_path. Compmlies these extractions as a dictionary, which is passed to outputDistributedExtractions_country for export.
    
    Params to pass are:
    
    slices       : a list of three slice objects definng start,stop,step for lat,long,month respectively.: e.g [slice(None,None,None), slice(None,None,None), slice(0,12,None)]
    a_lo,a_hi    : lower and upper age to predict for
    n_per        : how many realisations of the nugget are we simulating
    FileStartRel : number of first realisation present in the hdf5 file (in filename)
    FileEndRel   : number of last realisation (up to but not including) present in the hdf5 file (in filename)
    startRel     : number of first realisation WITHIN THIS FILE that we want to extract over (if ommited will start from 0)
    endRel       : number of last realisation WITHIN THIS FILE that we want to extract over (if ommited will use last realisation in file)
    ''' 

    ####TEMP
    #XXXa=r.Sys_time()
    #xxx1=xxx2=xxx3=xxx4=xxx5=xxx6=xxx7=xxx8=xxx9=xxx10=xxx11=xxx12=xxx13=xxx14=xxx15=xxx16=xxx17=xxx18=xxx19=xxx20=0
    ########

    # construct filepath for this realisation block, and define link
    filename = realisations_path
    filename = filename.replace('FILESTARTREL',str(FileStartRel))
    filename = filename.replace('FILEENDREL',str(FileEndRel))
    #checkAndBuildPaths(filename,VERBOSE=True,BUILD=False)
    hf = tb.openFile(filename)    
    hr = hf.root

    # define default start and end realisations WITHIN AND RELATIVE TO THIS FILE
    if startRel is None: startRel = 0 
    if endRel is None: endRel = hr.realizations.shape[0]
    
    # if either startRel or endRel are specified, run check that the hdf5 file contains sufficient realisations
    if ((startRel is None) & (endRel is None))==False:
        if((endRel - startRel)>hr.realizations.shape[0]):
            print 'ERROR!!! asking for '+str(endRel - startRel)+' realisations from block '+str(filename)+' that has only '+str(hr.realizations.shape[0])+' : EXITING!!!'
            return(-9999)

    #xxx1a = r.Sys_time()
    # define basic parameters
    slices = tuple(slices)     
    n_realizations = (endRel - startRel)
    n_rows=len(hr.lat_axis)
    n_cols=len(hr.lon_axis)
    N_facs = int(1e5)
    N_years = (slices[2].stop - slices[2].start)/12
    
    # define start and end rows for each iteation of loop (taking into account variable width of last remaining chunk)
    rowList=np.arange(0,n_rows)
    startRows=rowList[0:n_rows:rowsPerChunk]
    endRows = startRows+rowsPerChunk
    if endRows[-1]>(n_rows+1):
        endRows[-1]=(n_rows+1)
    NrowChunks = len(startRows)    

    # Get nugget variance and age-correction factors    
    V = hr.PyMCsamples.col('V')[:]    
    facs = mbgw.correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)    

    # open link to salb grid (masked to stable areas only) and population grid    
    salblim1km = tb.openFile(salblim1km_path, mode = "r")    
    grump1km = tb.openFile(grump1km_path, mode = "r")    

    # perform check that the number of rows and columns is the same in both 1km grids
    if len(salblim1km.root.lat) != len(grump1km.root.lat):
        print 'WARNING!! 1km row numbers do not correspond: salblim1km has '+str(len(salblim1km.root.lat))+' and grump1km has '+str(len(grump1km.root.lat))
    if len(salblim1km.root.long) != len(grump1km.root.long):
        print 'WARNING!! col numbers do not correspond: salblim1km has '+str(len(salblim1km.root.long))+' and grump1km has '+str(len(grump1km.root.long))

    # perform check that the number of rows and columns is in the correct ratio to those of input 5km grid
    if len(salblim1km.root.lat) != HiResLowResRatio*len(hr.lat_axis):
        print 'WARNING!! 1km and 5km row numbers do not correspond: salblim1km has '+str(len(salblim1km.root.lat))+' and 5km rows * HiResLowResRatio is '+str(HiResLowResRatio*len(hr.lat_axis))
    if len(salblim1km.root.long) != HiResLowResRatio*len(hr.lon_axis):
        print 'WARNING!! 1km and 5km col numbers do not correspond: salblim1km has '+str(len(salblim1km.root.long))+' and 5km cols * HiResLowResRatio is '+str(HiResLowResRatio*len(hr.lon_axis))

    # get list of unique salb IDs and count of pixels in each.. 
    # ..first check that Salb grid has been pre-examined using examineSalb and lists of unique IDs and N pixels exist, if not then re-run examineSalb
    try:
        uniqueSalb=fromfile(uniqueSalb_path,sep=",")
        pixelN=fromfile(pixelN_path,sep=",")
    except IOError:
        print 'WARNING!! files '+pixelN_path+" or "+uniqueSalb_path+" not found: running examineSalb"
        temp=examineSalb (salblim1km_path,ignore=np.array([-9999]))
        uniqueSalb=temp['uniqueSalb']
        pixelN=temp['count'] 

    #uniqueSalb=fromfile(uniqueSalb_path,sep=",")     
    #pixelN=fromfile(pixelN_path,sep=",")
    Nsalb=len(uniqueSalb)    

    # intialise empty arrays (e.g. 87 countries * N realisations) for mean PR in each country 
    countryMeanPRrel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per)     
    #countryBURDENrel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) 

    # intialise empty arrays (e.g. 87 countries * N realisations) for PAR and BURDEN in each class of each scheme in each country..housed in PAR dictionary PARdict and BURDEN dictionary BURDENdict
    Nschemes=len(breaksDict)    
    schemeNames=breaksDict.keys()    
    PARdict=cp.deepcopy(breaksDict)
    BURDENdict=cp.deepcopy(breaksDict)
    #xxx1 = xxx1 + (r.Sys_time() - xxx1a)
    
    #xxx2a = r.Sys_time()
    ## ..loop through each classification scheme 
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = PARdict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames) 

        # define additional sub-dictionaries to add to PARdict and BURDENdict to house arrays for PAR and BURDEN per class per scheme per country
        PAR = {}
        BURDEN = {}

        # .. for each class within each scheme..
        for cc in xrange(0,Nclasses):
            thisbreakname = breaknames[cc]
            
            # define two empty arrays for this scheme-class to house PAR and BURDEN per country realisations for each class
            blankarray_PAR = {thisbreakname: repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) }
            blankarray_BURDEN = {thisbreakname: repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) }

            # add these blank arrays to interim PAR and BURDEN dictionaries
            PAR.update(blankarray_PAR)
            BURDEN.update(blankarray_BURDEN)

        # add these sub-dictionary to PARdict for this scheme
        PAR = {'PAR':PAR}
        BURDEN = {'BURDEN':BURDEN}
        
        PARdict[scheme].update(PAR)
        BURDENdict[scheme].update(BURDEN)
    #xxx2 = xxx2 + (r.Sys_time() - xxx2a)

    # define a function object for later estimation of burden, basedon this grump1km row (after cnvertig to a vector)
    #ind1km = np.where(grump1km_ROW!=-99999999)
    #POPsurfaceVECTOR=grump1km_ROW[ind1km]
    BurdenPredictorObj = BurdenPredictor(hf_name=burdentrace_path, nyr=N_years, burn=0)

    # loop through each realisation
    for ii in xrange(0,n_realizations): #1:500 realisations n_realizations   
    
        # define which realisatin this relates to in global set from MCMC
        MCMCrel = startRel+ii 

        print "realisation "+str(MCMCrel)+" of set "+str(startRel)+" to "+str(endRel)+" ("+str(n_realizations)+" realisations)"    

        #xxx3a = r.Sys_time() 
        # Pull out parasite rate chunk (i.e. import n months of block)    
        tot_slice = (slice(MCMCrel,MCMCrel+1,None),) + slices    
        #f_chunk = hr.realizations[tot_slice][::-1,:,:,:].T   # hr.realizations=[rel,row,col,month]   #f_chunk = [month,col,row,rel]
        #f_chunk = f_chunk[:,:,:,0]                       #f_chunk = [month,col,row]

        # because pyTables is not working properly, manually loop through each month we want to extract and make a ST 3d matrix
        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = zeros(1*n_cols*n_rows*n_months).reshape(1,n_cols,n_rows,n_months)
        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,mm] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
        f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0]   

        ########TEMP###########
        #set missing vlaues in f block to 0
        #from scipy.io import write_array
        #write_array('/home/pwg/MBGWorld/extraction/temp_PRrel1.txt', f_chunk[0,:,:])
        print('NaNs in f_chunk:' +str(sum(isnan(f_chunk))))
        print('0s in f_chunk:' +str(sum(f_chunk==0)))
        f_chunk[isnan(f_chunk)]=0
        print('NaNs in f_chunk:' +str(sum(isnan(f_chunk))))
        print('0s in f_chunk:' +str(sum(f_chunk==0)))
        ####################################
        #xxx3 = xxx3 + (r.Sys_time() - xxx3a)

        #xxx4a = r.Sys_time() 
        # run check that there are no missing values in this f chunk
        if sum(isnan(f_chunk))>0:
            print "WARNING!! found "+str(sum(isnan(f_chunk)))+" NaN's in realisation "+str(MCMCrel)+" EXITING!!!"
            return(-9999)

        ## initialise arrays to house running mean PR whilst we loop through chunks and nugget draws..
        countryMeanPRrel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)
        #countryBURDENrel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)

        # initialise arrays for running PAR and running total burden whilst we loop through chunks and nugget draws  (housed in PARdic and BURDENdict, so simply ensure reset to zero for this realisation)..
        PARdict_ChunkRunning=cp.deepcopy(breaksDict)
        BURDENdict_ChunkRunning=cp.deepcopy(breaksDict)

        # ..loop through each classification scheme.. 
        for ss in xrange(0,Nschemes):
            scheme = schemeNames[ss]    
            breaknames = breaksDict[scheme]['BREAKNAMES']
            Nclasses=len(breaknames) 

            # define additional sub-dictionaries to add to PARdict_ChunkRunning to house temporary arrays
            PAR = {}
            BURDEN = {}

            # .. for each  class within each scheme..
            for cc in xrange(0,Nclasses):
                thisbreakname = breaknames[cc]
            
                # define an empty array for this scheme-class to house running PAR sum
                blankarray_PAR = {thisbreakname: repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per) }
                blankarray_BURDEN = {thisbreakname: repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per) }

                # add these blank arrays to interim PAR,MeanPR, and BURDEN dictionaries
                PAR.update(blankarray_PAR)
                BURDEN.update(blankarray_BURDEN)

            # add these sub-dictionaries to ChunkRunning dictionaries for this scheme
            PAR = {'PAR':PAR}
            BURDEN = {'BURDEN':BURDEN}
            PARdict_ChunkRunning[scheme].update(PAR)
            BURDENdict_ChunkRunning[scheme].update(BURDEN)
        #xxx4 = xxx4 + (r.Sys_time() - xxx4a)
        
        # loop through each row (or multiple rows in a slice) of 5km realisation grid..
        #n_slices = n_rows/rowsInslice5km

        timea = time.time()
        interimCnt=0 
        #for jj in xrange(0,n_rows): 
        for jj in xrange(0,NrowChunks): 
        
            # which rows of the 5km PR block are we dealing with in this iteration
            startRow=startRows[jj]
            endRow=endRows[jj]

            interimCnt=interimCnt+1
            if interimCnt==100:
                print('    on slice '+str(jj)+' of '+str(n_rows))
                print "slice time: "+str(time.time()-timea)
                timea=time.time()
                interimCnt=0
                    
            #xxx5a = r.Sys_time() 
            # get row of 5km PR surface accross all months in chunk  (assumes f_chunk is correct way up i.e. map view)
            #f_chunk_ROW = f_chunk[:,jj,:]
            f_chunk_ROW = f_chunk[:,startRow:endRow:1,:]
            
            # get corresponding 5 rows of 1km Salb and population surface (assumes they are correct way up i.e. map view)
            startRow1km=startRow*HiResLowResRatio
            endRow1km=endRow*HiResLowResRatio
            salblim1km_ROW = salblim1km.root.data[slice(startRow1km,endRow1km,1),:]
            grump1km_ROW = grump1km.root.data[slice(startRow1km,endRow1km,1),:]
            
            # define a blank array of zeroes of same size as 1km chunk - that will be duplicated for various uses later
            zeroChunk = zeros(np.product(grump1km_ROW.shape)).reshape(grump1km_ROW.shape)
            #xxx5 = xxx5 + (r.Sys_time() - xxx5a) 

            #plotMapPY(salblim1km.root.data[:,:],NODATA=-9999)
            #plotMapPY(salblim1km.root.data[slice(0,100,1),:],NODATA=-9999)
            #plotMapPY(f_chunk[0,:,:])
            #plotMapPY(f_chunk[0,slice(0,100,1),:])

            #xxx6a = r.Sys_time() 
            # how many unique salb IDs in these rows (after removing -9999 cells from sea and/or non-stable areas)
            uniqueSalb_ROW = unique(salblim1km_ROW)
            uniqueSalb_ROW = uniqueSalb_ROW[(uniqueSalb_ROW!=-9999)]
            Nsalb_ROW =len(uniqueSalb_ROW)

            # if we only have -9999 cells in this chunk, then can ignore and go to next chunk (ie. onto next jj)
            if Nsalb_ROW==0:
                continue

            # which rows on the all-country arrays correpsond to these salb IDs?
            salbLUT_ROW=array(r.match(uniqueSalb_ROW,uniqueSalb))
            #.. take care of special case of only 1 unique country ID in this slice - must ensure salbLUT_ROW is an array even if onluy one element
            if Nsalb_ROW==1:
                 salbLUT_ROW.shape=1
            salbLUT_ROW = salbLUT_ROW-1  # converts from r to python naming convention   
            #xxx6 = xxx6 + (r.Sys_time() - xxx6a)

            # do a pre-loop through each country present in this chunk, and define an ID matrix showing which pixels belong to it (but excluding non-stable areas)
            sumCountryID=0
            countryIDdict={}
            #xxx7a = r.Sys_time()
            for rr in xrange(0,Nsalb_ROW): 
                countryIDmatrix = cp.deepcopy(zeroChunk)
                countryIDmatrix[salblim1km_ROW==uniqueSalb_ROW[rr]]=1
                sumCountryID = sumCountryID + np.sum(countryIDmatrix)
                tmpdict={str(uniqueSalb_ROW[rr]):cp.deepcopy(countryIDmatrix)}
                countryIDdict.update(tmpdict)
            sumCountryID = sumCountryID + np.sum(salblim1km_ROW==-9999)    
            #xxx7 = xxx7 + (r.Sys_time() - xxx7a)
            
            # run check that sum of total pixels in all countries in this block match expected total
            if sumCountryID != np.product(salblim1km_ROW.shape):
                print "WARNING!! sum of pixels in each country in chunk "+str(jj)+" is "+str(sumCountryID)+" != expected total ("+str(np.product(salblim1km_ROW.shape))+")"

            # loop through n_per draws of the nugget..
            for kk in xrange(0,n_per):

                
                # add nugget component, apply inverse logit, apply age-correction factor
                #xxx8a = r.Sys_time()
                chunk = f_chunk_ROW + np.random.normal(loc=0, scale=np.sqrt(V[MCMCrel]), size=f_chunk_ROW.shape)
                chunk = pm.invlogit(chunk.ravel())
                chunk *= facs[np.random.randint(N_facs, size=np.prod(chunk.shape))]
                chunk = chunk.reshape(f_chunk_ROW.shape).squeeze()
                #xxx8 = xxx8 + (r.Sys_time() - xxx8a)

                # aggregate through time to obtain spatial-only array for this nugget-realisation
                #xxx9a = r.Sys_time()
                chunkTMEAN = atleast_2d(np.mean(chunk,0))
            
                # make a mappng vector for later conversion of arays of this dimension to a vector, and back again
                #ind5km = np.where(chunkTMEAN!=-99999999) 
                #xxx9 = xxx9 + (r.Sys_time() - xxx9a)

                # run check that this time-aggregated chunk has same spatial dimensions as time block
                #test=chunk.shape[1:]==chunkTMEAN.shape[0:]
                #if test==False:
                #    print("WARNING !!!!: spatial dimensions of time-aggregated block 'chunkTMEAN' do not match pre-aggregation 'chunk': EXITING!!")
                #    print ('chunk.shape[1:]: '+str(chunk.shape[1:]))
                #    print ('chunkTMEAN.shape[1:]: '+str(chunkTMEAN.shape[1:]))
                #    return(-9999)

                #xxx10a = r.Sys_time()
                # now expand the 5km PR chunk to match underlying 1km grid
                chunkExp = expandGridResPY(chunkTMEAN,HiResLowResRatio)
                #xxx10 = xxx10 + (r.Sys_time() - xxx10a)                

                # run check that this expanded block has correct dimensions
                test=chunkExp.shape==salblim1km_ROW.shape           
                if test==False:
                    print("WARNING !!!!: spatial dimensions of expanded 5km 'chunkExp' do not match 1km covariate chunk 'salblim1km_ROW': EXITING!! ")            
                    return(-9999)

                # run check that there are no PR==-9999 pixels (assigned to non-stable pixels in CS code) in stable areas on salblim1km
                testmatrix = cp.deepcopy(chunkExp)
                testmatrix[salblim1km_ROW == -9999] = 0
                if (np.sum(testmatrix == -9999) > 0):
                    print ("WARNING!!: ("+str(np.sum(testmatrix== -9999))+") null PR pixels (-9999) found in stable areas in rel "+str(ii)+" , row "+str(jj) )+ ": EXITING!!"
                    return(-9999)

                # obtain a burden surface for this chunk as a function of population and PR
                ## convert PRsurface and POPsurface to vectors before passing, then back=convert afterwards
                #PRsurfaceVECTOR=chunkTMEAN[ind5km]
                #print 'shape(chunkTMEAN): '+str(shape(chunkTMEAN))
                #print 'shape(grump1km_ROW): '+str(shape(grump1km_ROW)) 
                burdenChunk = BurdenPredictorObj.pr5km_pop1km(pr=chunkTMEAN.squeeze(),pop=grump1km_ROW,pop_pr_res=HiResLowResRatio)
                #burdenChunk = BurdenPredictorObj(pr=chunkTMEAN,pop=grump1km_ROW,pop_pr_res=HiResLowResRatio)
                #burdenChunk = grump1km_ROW*2
                #print 'shape(burdenChunk): '+str(shape(burdenChunk))
                #print 'shape(chunkTMEAN): '+str(shape(chunkTMEAN))+'\n'
                #burdenChunk =cp.deepcopy(chunkExp) # simply provides a template in correct format (1km) to populate with burdenChunkVECTOR
                #burdenChunk[ind1km]=burdenChunkVECTOR
                
                # create an ID matrix for this chunk for each endemicity class in each scheme                
                #xxx11a = r.Sys_time()
                classIDdict = cp.deepcopy(breaksDict)
                for ss in xrange(0,Nschemes): 
                    scheme = schemeNames[ss]   
                    breaknames = classIDdict[scheme]['BREAKNAMES']
                    breaks = classIDdict[scheme]['BREAKS']
                    Nclasses=len(breaknames) 

                    # define additional sub-dictionaries to add to classIDdict to house classID arrays for this chunk
                    classIDarrays = {}

                    # .. for each  class within each scheme..
                    #sumClassID=0
                    for cc in xrange(0,Nclasses):
                        thisbreakname=breaknames[cc]
                        
                        # define an ID matrix to identify those pixels in this chunk in in this class
                        classID = cp.deepcopy(zeroChunk)
                        classID[(chunkExp>=breaks[cc]) & (chunkExp < breaks[cc+1])]=1 
                        #sumClassID=sumClassID+sum(classID)
                        blankarray = {thisbreakname: cp.deepcopy(classID) }

                        # add these blank array to interim PAR dictionary
                        classIDarrays.update(blankarray)

                    # add these sub-dictionaries to PARdict for this scheme
                    classIDarrays = {'classIDarrays':classIDarrays}
                    classIDdict[scheme].update(classIDarrays)

                    # run check that sum of total pixels in all classes in this block match expected total
                    #if sumClassID != np.product(chunkExp.shape):
                    #    print "WARNING!! sum of pixels in each class in chunk "+str(jj)+" is "+str(sumClassID)+" != expected total ("+str(np.product(chunkExp.shape))+")"
                #xxx11 = xxx11 + (r.Sys_time() - xxx11a)

                # loop through each unique country in this chunk: calculate running mean PR, running total burden and PAR in each endemicity class in each scheme
                for rr in xrange(0,Nsalb_ROW):
                
                    thiscountry_salbLUT = salbLUT_ROW[rr] 

                    # get this country's ID matrix from countryIDdict dictionary
                    countryID = countryIDdict[str(uniqueSalb_ROW[rr])]

                    #xxx12a = r.Sys_time()
                    # calculate sum of PR in this country,convert to running mean using known country pixel count, and add to the relevant part of countryMeanPRrel_ChunkRunning
                    PRsum = np.sum(chunkExp*countryID)
                    countryMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk] = countryMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk]+(PRsum/pixelN[thiscountry_salbLUT])

                    #xxx13a = r.Sys_time()
                    # loop through schemes and classes to extract PAR for this country
                    for ss in xrange(0,Nschemes):
                        scheme = schemeNames[ss]    
                        breaknames = classIDdict[scheme]['BREAKNAMES']
                        breaks = classIDdict[scheme]['BREAKS']
                        Nclasses=len(breaknames) 
                        for cc in xrange(0,Nclasses):
                            thisbreakname = breaknames[cc] 
                            # get this classes's ID matrix from classIDdict dictionary
                            classID = classIDdict[scheme]['classIDarrays'][thisbreakname]

                            # define ID matrix for this country AND this class
                            countryClassID = countryID*classID    
                            
                            # calculate sum of population in this class and country and add this sum to the relevant part of PARdict_ChunkRunning
                            PARtemp=grump1km_ROW * countryClassID
                            PARsum = np.sum(PARtemp)
                            PARdict_ChunkRunning[scheme]['PAR'][thisbreakname][thiscountry_salbLUT,kk] = PARdict_ChunkRunning[scheme]['PAR'][thisbreakname][thiscountry_salbLUT,kk] + PARsum

                            # similarly, calculate sum of burden in this country, and add to relevant part of countryBURDENrel_ChunkRunning 
                            BURDENtemp = burdenChunk*countryClassID
                            BURDENsum = np.sum(BURDENtemp)
                            BURDENdict_ChunkRunning[scheme]['BURDEN'][thisbreakname][thiscountry_salbLUT,kk] = BURDENdict_ChunkRunning[scheme]['BURDEN'][thisbreakname][thiscountry_salbLUT,kk] + BURDENsum
                            #xxx12 = xxx12 + (r.Sys_time() - xxx12a)                            

                            # run check for nonsensical negative PAR or BURDEN value
                            if PARsum<0.:
                                print('WARNING!! Negative PAR found - check input population grid. EXITING')
                                return(-9999)
                                
                            if BURDENsum<0.:
                                print('WARNING!! Negative BURDEN found - check input population grid. EXITING')
                                return(-9999)
                    #xxx13 = xxx13 + (r.Sys_time() - xxx13a)

        # copy mean PR values for these 500 nugget draws to the main arrays housing all realisations
        #xxx17a = r.Sys_time()
        countryMeanPRrel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryMeanPRrel_ChunkRunning
        #countryBURDENrel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryBURDENrel_ChunkRunning

        # loop through class schemes
        for ss in xrange(0,Nschemes):
            scheme =  schemeNames[ss]   
            breaknames = breaksDict[scheme]['BREAKNAMES']
            Nclasses=len(breaknames) 

            # .. for each class within each scheme..
            for cc in xrange(0,Nclasses):
                thisbreakname = breaknames[cc]
                
                # copy PAR values for these 500 nugget draws to main arrays housing all realisations
                PARdict[scheme]['PAR'][thisbreakname][:,slice(ii*n_per,(ii*n_per)+n_per,1)] = PARdict_ChunkRunning[scheme]['PAR'][thisbreakname] 

                # copy BURDEN values for these 500 nugget draws to main arrays housing all realisations
                BURDENdict[scheme]['BURDEN'][thisbreakname][:,slice(ii*n_per,(ii*n_per)+n_per,1)] = BURDENdict_ChunkRunning[scheme]['BURDEN'][thisbreakname] 
        #xxx17 = xxx17 + (r.Sys_time() - xxx17a)
        
    #####TEMP
    #XXX=r.Sys_time()-XXXa
    #print 'TOTAL time         : '+str(XXX)+' ('+str((XXX/XXX)*100)+'%)'
    #print '1   : '+str(xxx1)+' ('+str((xxx1/XXX)*100)+'%)'
    #print '2   : '+str(xxx2)+' ('+str((xxx2/XXX)*100)+'%)'
    #print '3   : '+str(xxx3)+' ('+str((xxx3/XXX)*100)+'%)'
    #print '4   : '+str(xxx4)+' ('+str((xxx4/XXX)*100)+'%)'
    #print 'assigning fchunk   : '+str(xxx5)+' ('+str((xxx5/XXX)*100)+'%)'
    #print '6   : '+str(xxx6)+' ('+str((xxx6/XXX)*100)+'%)'
    #print '7   : '+str(xxx7)+' ('+str((xxx7/XXX)*100)+'%)'
    #print 'fchunk to PR       : '+str(xxx8)+' ('+str((xxx8/XXX)*100)+'%)'
    #print 'time aggregation   : '+str(xxx9)+' ('+str((xxx9/XXX)*100)+'%)'
    #print '10   : '+str(xxx10)+' ('+str((xxx10/XXX)*100)+'%)'
    #print '11   : '+str(xxx11)+' ('+str((xxx11/XXX)*100)+'%)'
    #print '12   : '+str(xxx12)+' ('+str((xxx12/XXX)*100)+'%)'
    #print '13   : '+str(xxx13)+' ('+str((xxx13/XXX)*100)+'%)'
    #print '14   : '+str(xxx14)+' ('+str((xxx14/XXX)*100)+'%)'
    #print '15   : '+str(xxx15)+' ('+str((xxx15/XXX)*100)+'%)'
    #print '16   : '+str(xxx16)+' ('+str((xxx16/XXX)*100)+'%)'
    #print '17   : '+str(xxx17)+' ('+str((xxx17/XXX)*100)+'%)'


    # define absolute startRel and endRel for appropriate output file suffixes
    startRelOUT = FileStartRel+startRel
    endRelOUT = FileStartRel+endRel
    
    # define dictionary to pass to outputDistributedExtractions_country containing:
    # 1. Array of mean PR values per country per realisation (countryMeanPRrel)
    # 2. Array of burden totals per country per realisation (countryBURDENrel)
    # 3. dictionary of PAR values per country per realisation, along with classification scheme metadata (PARdict)

    returnDict={"startRel":startRelOUT,"endRel":endRelOUT,"countryMeanPRrel":countryMeanPRrel,"BURDENdict":BURDENdict,"PARdict":PARdict}
    
    # export extracted arrays 
    outputDistributedExtractions_country(returnDict)
    
    #return(returnDict)
#############################################################################################################################################
def outputDistributedExtractions_country(dict): 

    ''''
    Takes a dictionary which is output from extractSummaries_country and exports 
    .txt files of mean PR, burden, and PAR extractions by country.
    
    Params to pass are:
    
    dict      : output from outputDistributedExtractions_country 
    '''

    # check for error output from extractSummaries_country due to NaNs
    if dict == -9999:
        print "WARNING!! recieved error output from extractSummaries_country() - will not run outputDistributedExtractions_country()" 
        return(-9999)

    # define which realisations we are dealing with
    startRel=dict['startRel']
    endRel=dict['endRel']

    # construct file suffix indicating realisations in question
    relSuff = '_r'+str(startRel)+'to'+str(endRel)

    # export arrays of mean PR per country per realisation
    np.savetxt(exportPathDistributed_country+'meanPR_country'+relSuff+'.txt', dict['countryMeanPRrel'])
    #np.savetxt(exportPathDistributed_country+'BURDEN_country'+relSuff+'.txt', dict['countryBURDENrel'])

    Nschemes=len(breaksDict)    
    schemeNames=breaksDict.keys()    

    # loop through classification schemes and clases wihtin them and export array of PAR for each
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = breaksDict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames)

        for cc in xrange(0,Nclasses):

            # construct class and scheme suffix
            classSuff = '_'+schemeNames[ss]+'_'+breaknames[cc]

            # export PAR array for this scheme-class
            np.savetxt(exportPathDistributed_country+'PAR_country'+classSuff+relSuff+'.txt', dict['PARdict'][schemeNames[ss]]['PAR'][breaknames[cc]]) 

            # export BURDEN array for this scheme-class
            np.savetxt(exportPathDistributed_country+'BURDEN_country'+classSuff+relSuff+'.txt', dict['BURDENdict'][schemeNames[ss]]['BURDEN'][breaknames[cc]]) 
#############################################################################################################################################
def extractSummaries_perpixel (slices,a_lo,a_hi,n_per,FileStartRel,FileEndRel,totalN,startRel=None,endRel=None,BURDEN=False):

    '''
    Takes an hdf block of one or more realisations of f, and calculates running arrays that can subsequently be combined using combineDistribExtractions_perpixel()
    to calculate per-pixel summaries (mean and StdDev) of posterior distributions of PR and burden. Also calculates posterior
    probability of membership to every class in every scheme (specified by breaksDict), most likely class in each scheme, and
    probability of membership to that most likely class. Export gzipped arrays for subsequent import and combining in combineDistribExtractions_perpixel() 
        
    Params to pass are:
    
    slices       : a list of three slice objects defining start,stop,step for lat,long,month respectively.: e.g [slice(None,None,None), slice(None,None,None), slice(0,12,None)]
    a_lo,a_hi    : lower and upper age to predict for
    n_per        : how many realisations of the nugget are we simulating
    FileStartRel : number of first realisation present in the hdf5 file (in filename)
    FileEndRel   : number of last realisation (up to but not including) present in the hdf5 file (in filename)
    totalN       : the total number of realisations (i.e. denominator in posterior mean) i.e n_per * n_realizations
    startRel     : number of first realisation WITHIN THIS FILE that we want to extract over (if ommited will start from 0)
    endRel       : number of last realisation WITHIN THIS FILE that we want to extract over (if ommited will use last realisation in file)
    BURDEN       : do we want to perform calculations for burden - default is no.
    ''' 

    # construct filepath for this realisation block, and define link
    filename = realisations_path
    filename = filename.replace('FILESTARTREL',str(FileStartRel))
    filename = filename.replace('FILEENDREL',str(FileEndRel))
    hf = tb.openFile(filename)    
    hr = hf.root

    # define default start and end realisations WITHIN AND RELATIVE TO THIS FILE
    if startRel is None: startRel = 0 
    if endRel is None: endRel = hr.realizations.shape[0]
    
    # if either startRel or endRel are specified, run check that the hdf5 file contains sufficient realisations
    if ((startRel is None) & (endRel is None))==False:
        if((endRel - startRel)>hr.realizations.shape[0]):
            print 'ERROR!!! asking for '+str(endRel - startRel)+' realisations from block '+str(filename)+' that has only '+str(hr.realizations.shape[0])+' : EXITING!!!'
            return(-9999)

    #xxx1a = r.Sys_time()
    # define basic parameters
    slices = tuple(slices)     
    n_realizations = (endRel - startRel)
    n_rows=len(hr.lat_axis)
    n_cols=len(hr.lon_axis)
    N_facs = int(1e5)
    N_years = (slices[2].stop - slices[2].start)/12

    # Get nugget variance and age-correction factors    
    V = hr.PyMCsamples.col('V')[:]    
    facs = mbgw.correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)    

    # if we are extracting burden summaries, import 5km population grid
    if BURDEN==True:
        grump5km = tb.openFile(grump5km_path)
        
        # define a function object for later estimation of burden, basedon this grump1km row (after cnvertig to a vector)
        #ind = np.where(grump5km.root.data[:,:]!=-99999999)
        #POPsurfaceVECTOR=grump5km.root.data[:,:][ind]
        BurdenPredictorObj = BurdenPredictor(hf_name=burdentrace_path, nyr=N_years, burn=0) 

    # define a blank array of zeroes of same size as a single monthly map - that will be duplicated for various uses later
    zeroMap = np.zeros(n_rows*n_cols).reshape(n_rows,n_cols)

    # initialise zero matrices that will house running totals
    meanPR = cp.deepcopy(zeroMap)
    meanPR2 = cp.deepcopy(zeroMap)
    if BURDEN==True:
        meanBUR = cp.deepcopy(zeroMap)    
        meanBUR2 = cp.deepcopy(zeroMap)

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

        print "realisation "+str(FileStartRel+MCMCrel)+" of set "+str(FileStartRel)+" to "+str(FileEndRel)+" ("+str(n_realizations)+" realisations)"    

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
        #if sum(sum(sum(np.isnan(f_chunk))))>0:
        if np.isnan(f_chunk).any()==True:
        
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
            #exportAscii(chunkTMEAN,exportPathDistributed_perpixel+"chunkTmean.asc",hdrDict)
            #return()
            
            # increment runing mean PR matrices 
            meanPR = meanPR + (chunkTMEAN/totalN)
            meanPR2 = meanPR2 + (np.square(chunkTMEAN)/totalN)
            
            # get burden realisation for this PR and increment running burden matrix
            if BURDEN==True:

                ## convert PRsurface to vector before passing, then back=convert afterwards
                #PRsurfaceVECTOR=chunkTMEAN[ind]
                burdenChunk = BurdenPredictorObj.pr5km_pop5km(pr=chunkTMEAN,pop=grump5km.root.data[:,:])
                #burdenChunk =cp.deepcopy(chunkTMEAN) # simply provides a template in correct format to populate with burdenChunkVECTOR
                #burdenChunk[ind]=burdenChunkVECTOR

                meanBUR = meanBUR + (burdenChunk/totalN)
                meanBUR2 = meanBUR2 + (np.square(burdenChunk)/totalN)
            
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

                    # update class membership running probability array
                    PCMdict[scheme]['PCM'][thisbreakname] = PCMdict[scheme]['PCM'][thisbreakname] + (classID/totalN)

            
    # export running arrays for this set of realisations

    ## construct file suffix indicating realisations in question
    startRelOUT = FileStartRel+startRel
    endRelOUT = FileStartRel+endRel
    relSuff = '_r'+str(startRelOUT)+'to'+str(endRelOUT)

    ## export running meanPR and meanPR2 array
    np.savetxt(exportPathDistributed_perpixel+"meanPR_perpixel"+relSuff+".gz",meanPR)
    np.savetxt(exportPathDistributed_perpixel+"meanPR2_perpixel"+relSuff+".gz",meanPR2)


    if BURDEN==True:
        ## export running meanBUR and meanBUR2 array
        np.savetxt(exportPathDistributed_perpixel+"meanBUR_perpixel"+relSuff+".gz",meanBUR)
        np.savetxt(exportPathDistributed_perpixel+"meanBUR2_perpixel"+relSuff+".gz",meanBUR2)

    ## for each classification scheme, export running PCM arrays for each scheme/class
    for ss in xrange(0,Nschemes):             
        scheme=schemeNames[ss]                
        breaknames = PCMdict[scheme]['BREAKNAMES']             
        Nclasses=len(breaknames)

        # .. for each class within each scheme..            
        for cc in xrange(0,Nclasses):            
            thisbreakname = breaknames[cc]

            # whilst at this loop location, export PCM for this scheme/class as ascii
            np.savetxt(exportPathDistributed_perpixel+'PCM_perpixel_'+scheme+'_'+thisbreakname+relSuff+".gz",PCMdict[scheme]['PCM'][thisbreakname])

    return()  
#############################################################################################################################################
