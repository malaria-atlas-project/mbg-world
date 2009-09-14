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
#from mbgw.joint_simulation import *
from math import sqrt
from mbgw import correction_factors
from map_utils import getAsciiheaderFromTemplateHDF5
from map_utils import exportAscii
from map_utils import checkAndBuildPaths
#from IPython.Debugger import Pdb
#pdb = Pdb(color_scheme='Linux')

rowsPerChunk=1

# import R function
r.source('extract_Rlib.R')
expandGridResPY=r['expandGridRes']

r.source('/home/pwg/rc/PRutils.R')

# import parameters from param file
from extract_params import *
 
#############################################################################################################################################
def examineSalb (salb, grump=None,salbpop_path=None, pixarea=None,salbarea_path=None,uniqueSalb_path=None,pixelN_path=None,ignore=None):

    ''''
    Takes an input raster with integers defining unique spatial areas (e.g. the salb grid)
    and calculates two vectors: uniqueSalb is a list of unique integer IDs present in this raster,
    and count is the corresponding number of pixels in each. These are returned as a dictionary
    and are optionally exported to txt files.
    
    Params to pass are:
    
    salb            : either an hdf5 filepath or a numpy array object containing a raster of unique spatial unit IDs
    grump           : either an hdf5 filepath or a numpy array object containing a population grid. If present, will calculate pop per salb unit
    pixarea         : either an hdf5 filepath or a numpy array object containing a pixel area grid. If present, will calculate area per salb unit
    uniqueSalb_path : filepath to export uniqueSalb - no export if ommited
    pixelN_path     : filepath to export pixel count - no export if ommited
    salbpop_path    : filepath to export population count per salb unit - no export if ommited
    salbarea_path   : filepath to export area per salb unit - no export if ommited
    ignore          : a list containing values of salb to ignore from output lists (e.g. -9999)
    '''

    # if the salb file is specified with a filepath, import it  
    if type(salb) == str:
        salb = tb.openFile(salb, mode = "r")
    
    # if either the grump or pix area grids are sepcified with a filepath, import them, and check shapes match salb
    if type(grump) == str:
        grump = tb.openFile(grump, mode = "r")
        if (shape(salb.root.data)!=shape(grump.root.data)): raise RuntimeError ('Shape of salb grid '+str(shape(salb.root.data))+' != shape of grump grid '+str(shape(grump.root.data)))
    if type(pixarea) == str:
        pixarea = tb.openFile(pixarea, mode = "r")
        if (shape(pixarea.root.data)!=shape(grump.root.data)): raise RuntimeError ('Shape of salb grid '+str(shape(salb.root.data))+' != shape of pixarea grid '+str(shape(pixarea.root.data)))

#    from rpy import *
#    r.source('~/mbg-world/mbgw-scripts/extract_Rlib.R')
#    plotMapPY=r['plotMap']
#    r.par(mfrow=(1,3))
#    plotMapPY(salb.root.data[:,:],NODATA=-9999)
#    plotMapPY(grump.root.data[:,:],NODATA=-9999)
#    plotMapPY(pixarea.root.data[:,:],NODATA=-9999)
        
    # initialise some values
    nrow=shape(salb.root.data)[0]
    uniqueSalb=[-9999]
    count=[0]
    if (grump is not None): pop=[0]
    if (pixarea is not None): area=[0]    

    for i in xrange(0,nrow):
        #print(r.paste("row",i,"of",nrow))
        salb_chunk=salb.root.data[i,:]
        uniqueSalbCHUNK=unique(salb_chunk)
        NuniqueSalbCHUNK=len(uniqueSalbCHUNK)
        countCHUNK = zeros(NuniqueSalbCHUNK)
        if (grump is not None):
            popCHUNK = zeros(NuniqueSalbCHUNK)
            grump_chunk = grump.root.data[i,:]
        if (pixarea is not None):
            areaCHUNK = zeros(NuniqueSalbCHUNK)
            pixarea_chunk = pixarea.root.data[i,:]
            
        for j in xrange(0,NuniqueSalbCHUNK):
           
            where_salb = where(salb_chunk==uniqueSalbCHUNK[j])
            countCHUNK[j]=np.sum(salb_chunk==uniqueSalbCHUNK[j])

            if uniqueSalbCHUNK[j] in uniqueSalb:
                index=where(uniqueSalb==uniqueSalbCHUNK[j])[0]
                count[index]=count[index]+countCHUNK[j]
                if grump is not None:
                    popCHUNK[j] = np.sum(grump_chunk[where_salb])
                    pop[index]=pop[index]+popCHUNK[j]
                if pixarea is not None:
                    areaCHUNK[j] = np.sum(pixarea_chunk[where_salb])
                    area[index]=area[index]+areaCHUNK[j]
            else:
                uniqueSalb=append(uniqueSalb,uniqueSalbCHUNK[j]) 
                count=append(count,countCHUNK[j])
                if grump is not None:pop=append(pop,popCHUNK[j])
                if pixarea is not None:area=append(area,areaCHUNK[j])

    # optionally remove stated values (e.g 0 or -9999) from list of unique salbs
    if ignore is not None:
        for i in xrange(0,len(ignore)):
            count = count[uniqueSalb!=ignore[i]] 
            if grump is not None:pop = pop[uniqueSalb!=ignore[i]] 
            if pixarea is not None:area = area[uniqueSalb!=ignore[i]] 
            uniqueSalb = uniqueSalb[uniqueSalb!=ignore[i]]

    # optionally export to file
    if uniqueSalb_path is not None:
        uniqueSalb.tofile(uniqueSalb_path,sep=",")
    if pixelN_path is not None:    
        count.tofile(pixelN_path,sep=",")
    if salbpop_path is not None:    
        pop.tofile(salbpop_path,sep=",")
    if salbarea_path is not None:    
        area.tofile(salbarea_path,sep=",")

#    from rpy import *
#    r.par(mfrow=(1,2))
#    r.plot(count,area,xlab="count",ylab="area",main="")
#    r.plot(count,pop,xlab="count",ylab="pop",main="")
                
    # return unique list and corresponding pixel count as a dictionary
    returnDict={'uniqueSalb':uniqueSalb,'count':count}
    if grump is not None: returnDict.update({'pop':pop})
    if pixarea is not None: returnDict.update({'area':area})
    return returnDict
#############################################################################################################################################
def extractSummaries_country(slices,a_lo,a_hi,n_per,FileStartRel,FileEndRel,startRel=None,endRel=None, do_AREALMEANPR=False, do_POPMEANPR=False, do_PAR=False,do_BURDEN=False, do_AREALMEANRo=False, do_POPMEANRo=False):

    '''
    Takes an hdf block of one or more realisation of f, and calculates mean PR, total burden, and population at risk (PAR) in each
    unique spatial unit (e.g. country) specified in the salblim_path file. Also requires population surface specified by 
    grump_path. Compmlies these extractions as a dictionary, which is passed to outputDistributedExtractions_country for export.
    
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
    
    # define paths to input files according to specified resolution
    if (HiResLowResRatio_PERCOUNTRY==1):
        salblim_path=salblim5km_path
        salb_path=salb5km_path
        grump_path=grump5km_path
        pixarea_path=pixarea5km_path
        limbnry_path=lim5kmbnry_path
    if (HiResLowResRatio_PERCOUNTRY==5):
        salblim_path=salblim1km_path
        salb_path=salb1km_path
        grump_path=grump1km_path
        pixarea_path=pixarea1km_path
        limbnry_path=lim1kmbnry_path
    HiResLowResRatio=HiResLowResRatio_PERCOUNTRY
    
    # check we are actualy asking for one of the aggregated products
    if ((do_BURDEN is False) & (do_AREALMEANPR is False) & (do_POPMEANPR is False) & (do_POPMEANRo is False) & (do_AREALMEANRo is False) & (do_PAR is False)):
        print 'ERROR!!! Not asking for any aggregated products from extract_summaries_country, EXITING!!'
        return(-9999)

    # construct filepath for this realisation block, and define link
    filename = realisations_path
    filename = filename.replace('FILESTARTREL',str(FileStartRel))
    filename = filename.replace('FILEENDREL',str(FileEndRel))
    #checkAndBuildPaths(filename,VERBOSE=True,BUILD=False)
    hf = tb.openFile(filename)    
    hr = hf.root
    
    FileStartRel = int(FileStartRel)
    FileEndRel = int(FileEndRel)

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
    N_years = (slices[2].stop - slices[2].start)/12.
    
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

    # open link to salb grid (masked to stable areas only) and optionally population and pixarea grids    
    salblim = tb.openFile(salblim_path, mode = "r")     

    if (do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR):
        grump = tb.openFile(grump_path, mode = "r")    

        # perform check that the number of rows and columns is the same in both salb and grump grids
        if len(salblim.root.lat) != len(grump.root.lat):
            print 'WARNING!! row numbers do not correspond: salblim has '+str(len(salblim.root.lat))+' and grump has '+str(len(grump.root.lat))
        if len(salblim.root.long) != len(grump.root.long):
            print 'WARNING!! col numbers do not correspond: salblim has '+str(len(salblim.root.long))+' and grump has '+str(len(grump.root.long))

    if (do_AREALMEANPR | do_AREALMEANRo):
        pixarea = tb.openFile(pixarea_path, mode = "r")    

        # perform check that the number of rows and columns is the same in both salb and pixarea grids
        if len(salblim.root.lat) != len(pixarea.root.lat):
            print 'WARNING!! row numbers do not correspond: salblim has '+str(len(salblim.root.lat))+' and pixarea has '+str(len(pixarea.root.lat))
        if len(salblim.root.long) != len(pixarea.root.long):
            print 'WARNING!! col numbers do not correspond: salblim has '+str(len(salblim.root.long))+' and pixarea has '+str(len(pixarea.root.long))

    # perform check that the number of rows and columns is in the correct ratio to those of input 5km PR grid
    if len(salblim.root.lat) != HiResLowResRatio*len(hr.lat_axis):
        print 'WARNING!! PR grid and salb grid do not correspond: salblim has '+str(len(salblim.root.lat))+' and 5km PR grid rows * HiResLowResRatio is '+str(HiResLowResRatio*len(hr.lat_axis))
    if len(salblim.root.long) != HiResLowResRatio*len(hr.lon_axis):
        print 'WARNING!! PR grid and salb grid do not correspond: salblim has '+str(len(salblim.root.long))+' and 5km PR grid cols * HiResLowResRatio is '+str(HiResLowResRatio*len(hr.lon_axis))

    # get list of unique salb IDs and count of pixels in each.. 
    # ..first check that Salb grid has been pre-examined using examineSalb and lists of unique IDs and N pixels exist, if not then re-run examineSalb
    try:
        uniqueSalb=fromfile(uniqueSalb_path,sep=",")
        pixelN=fromfile(pixelN_path,sep=",")
        if (do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR): totalPop = fromfile(salbpop_path,sep=",")
        if (do_AREALMEANPR | do_AREALMEANRo): totalArea = fromfile(salbarea_path,sep=",")
    except IOError:
        print 'WARNING!! some necessary files from examineSalb were not found: re-running examineSalb'
        if (((do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR)==True) & ((do_AREALMEANPR | do_AREALMEANRo)==False)):
            temp = examineSalb(salb=salblim_path, grump=grump_path,salbpop_path=salbpop_path, pixarea=None,salbarea_path=None,uniqueSalb_path=uniqueSalb_path,pixelN_path=pixelN_path,ignore=np.array([-9999]))
            totalPop=temp['totalPop']
            
        if (((do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR)==False) & ((do_AREALMEANPR | do_AREALMEANRo)==True)):
            temp = examineSalb(salb=salblim_path, grump=None,salbpop_path=None, pixarea=pixarea_path,salbarea_path=salbarea_path,uniqueSalb_path=uniqueSalb_path,pixelN_path=pixelN_path,ignore=np.array([-9999]))
            totalArea=temp['totalArea']

        if (((do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR)==True) & ((do_AREALMEANPR | do_AREALMEANRo)==True)):
            temp = examineSalb(salb=salblim_path, grump=grump_path,salbpop_path=salbpop_path, pixarea=pixarea_path,salbarea_path=salbarea_path,uniqueSalb_path=uniqueSalb_path,pixelN_path=pixelN_path,ignore=np.array([-9999])) 
            totalPop=temp['totalPop']
            totalArea=temp['totalArea']
            
        uniqueSalb=temp['uniqueSalb']
        pixelN=temp['count'] 

    Nsalb=len(uniqueSalb)    

    # intialise empty arrays (e.g. 87 countries * N realisations) for mean PR in each country 
    if (do_AREALMEANPR): countryAreaMeanPRrel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per)     
    if (do_POPMEANPR): countryPopMeanPRrel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per)   
    if (do_AREALMEANRo): countryAreaMeanRorel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) 
    if (do_POPMEANRo): countryPopMeanRorel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per)   

    # intialise empty arrays (e.g. 87 countries * N realisations) for PAR and/or BURDEN in each class of each scheme in each country..housed in PAR dictionary PARdict and BURDEN dictionary BURDENdict
    if (do_PAR | do_BURDEN):
        Nschemes=len(breaksDict)    
        schemeNames=breaksDict.keys()    
        if (do_PAR): PARdict=cp.deepcopy(breaksDict)
        if (do_BURDEN): BURDENdict=cp.deepcopy(breaksDict)

        ## ..loop through each classification scheme 
        for ss in xrange(0,Nschemes): 
            scheme=schemeNames[ss]   
            breaknames = breaksDict[scheme]['BREAKNAMES']
            Nclasses=len(breaknames) 

            # define additional sub-dictionaries to add to PARdict and BURDENdict to house arrays for PAR and BURDEN per class per scheme per country
            if (do_PAR): PAR = {}
            if (do_BURDEN):BURDEN = {}

            # .. for each class within each scheme..
            for cc in xrange(0,Nclasses):
                thisbreakname = breaknames[cc]
                
                # define two empty arrays for this scheme-class to house PAR and BURDEN per country realisations for each class
                if (do_PAR): blankarray_PAR = {thisbreakname: repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) }
                if (do_BURDEN):blankarray_BURDEN = {thisbreakname: repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) }

                # add these blank arrays to interim PAR and BURDEN dictionaries
                if (do_PAR): PAR.update(blankarray_PAR)
                if (do_BURDEN):BURDEN.update(blankarray_BURDEN)

            # add these sub-dictionary to PARdict for this scheme
            if (do_PAR):PAR = {'PAR':PAR}
            if (do_BURDEN):BURDEN = {'BURDEN':BURDEN}
            
            if (do_PAR):PARdict[scheme].update(PAR)
            if (do_BURDEN):BURDENdict[scheme].update(BURDEN)

        # define a function object for later estimation of burden, based on this grump row (after cnvertig to a vector)
        if (do_BURDEN): BurdenPredictorObj = BurdenPredictor(hf_name=burdentrace_path, nyr=N_years, burn=0)

    # loop through each realisation
    for ii in xrange(0,n_realizations): #1:500 realisations n_realizations   
    
        # define which realisatin this relates to in global set from MCMC
        MCMCrel = startRel+ii 

        print "realisation "+str(MCMCrel)+" of set "+str(startRel)+" to "+str(endRel)+" ("+str(n_realizations)+" realisations)"    

        # Pull out parasite rate chunk (i.e. import n months of block)    
        tot_slice = (slice(MCMCrel,MCMCrel+1,None),) + slices    

        # because pyTables is not working properly, manually loop through each month we want to extract and make a ST 3d matrix
        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
        subsetmonth=0 
        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
            subsetmonth=subsetmonth+1
        f_chunk = f_chunk.squeeze()

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

        # run check that there are no missing values in this f chunk
        if sum(isnan(f_chunk))>0:
            raise RuntimeError ("Found "+str(sum(isnan(f_chunk)))+" NaN's in realisation "+str(MCMCrel)+" EXITING!!!")

        ## initialise arrays to house running mean PR whilst we loop through chunks and nugget draws..
        if (do_AREALMEANPR): countryAreaMeanPRrel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)
        if (do_POPMEANPR): countryPopMeanPRrel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)
        if (do_AREALMEANRo): countryAreaMeanRorel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)
        if (do_POPMEANRo): countryPopMeanRorel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)
        #countryBURDENrel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)

        if (do_PAR | do_BURDEN):

            # initialise arrays for running PAR and running total burden whilst we loop through chunks and nugget draws  (housed in PARdic and BURDENdict, so simply ensure reset to zero for this realisation)..
            if(do_PAR):PARdict_ChunkRunning=cp.deepcopy(breaksDict)
            if(do_BURDEN):BURDENdict_ChunkRunning=cp.deepcopy(breaksDict)

            # ..loop through each classification scheme.. 
            for ss in xrange(0,Nschemes):
                scheme = schemeNames[ss]    
                breaknames = breaksDict[scheme]['BREAKNAMES']
                Nclasses=len(breaknames) 

                # define additional sub-dictionaries to add to PARdict_ChunkRunning to house temporary arrays
                if(do_PAR):PAR = {}
                if(do_BURDEN):BURDEN = {}

                # .. for each  class within each scheme..
                for cc in xrange(0,Nclasses):
                    thisbreakname = breaknames[cc]
                
                    # define an empty array for this scheme-class to house running PAR sum
                    if(do_PAR):blankarray_PAR = {thisbreakname: repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per) }
                    if(do_BURDEN):blankarray_BURDEN = {thisbreakname: repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per) }

                    # add these blank arrays to interim PAR,MeanPR, and BURDEN dictionaries
                    PAR.update(blankarray_PAR)
                    if(do_BURDEN):BURDEN.update(blankarray_BURDEN)

                # add these sub-dictionaries to ChunkRunning dictionaries for this scheme
                if(do_PAR):PAR = {'PAR':PAR}
                if(do_BURDEN):BURDEN = {'BURDEN':BURDEN}
                if(do_PAR):PARdict_ChunkRunning[scheme].update(PAR)
                if(do_BURDEN):BURDENdict_ChunkRunning[scheme].update(BURDEN)
        
        # loop through each row (or multiple rows in a slice) of 5km realisation grid..
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
            f_chunk_ROW = f_chunk[startRow:endRow:1,:,:]
            
            # get corresponding 5km rows of 1km (or 5km) Salb and optionally population and pixarea surfaces (assumes they are correct way up i.e. map view)
            startRow1km=startRow*HiResLowResRatio  # NB '1km' suffix is retained here, to distinguish from rows of PR surface, even though the slab/grump etc surfaces can also be at 5km
            endRow1km=endRow*HiResLowResRatio
            salblim_ROW = salblim.root.data[slice(startRow1km,endRow1km,1),:]
            if(do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR): grump_ROW = grump.root.data[slice(startRow1km,endRow1km,1),:]
            if(do_AREALMEANPR | do_AREALMEANRo): pixarea_ROW = pixarea.root.data[slice(startRow1km,endRow1km,1),:]            
            
            # define a blank array of zeroes of same size as fine-res chunk - that will be duplicated for various uses later
            zeroChunk = zeros(np.product(salblim_ROW.shape)).reshape(salblim_ROW.shape)

            # how many unique salb IDs in these rows (after removing -9999 cells from sea and/or non-stable areas)
            uniqueSalb_ROW = unique(salblim_ROW)
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

            # do a pre-loop through each country present in this chunk, and define an ID matrix showing which pixels belong to it (but excluding non-stable areas)
            sumCountryID=0
            countryIDdict={}
            for rr in xrange(0,Nsalb_ROW): 
                countryIDmatrix = cp.deepcopy(zeroChunk)
                countryIDmatrix[salblim_ROW==uniqueSalb_ROW[rr]]=1
                sumCountryID = sumCountryID + np.sum(countryIDmatrix)
                tmpdict={str(uniqueSalb_ROW[rr]):cp.deepcopy(countryIDmatrix)}
                countryIDdict.update(tmpdict)
            sumCountryID = sumCountryID + np.sum(salblim_ROW==-9999)    
            
            # run check that sum of total pixels in all countries in this block match expected total
            if sumCountryID != np.product(salblim_ROW.shape):
                print "WARNING!! sum of pixels in each country in chunk "+str(jj)+" is "+str(sumCountryID)+" != expected total ("+str(np.product(salblim_ROW.shape))+")"

            # loop through n_per draws of the nugget..
            for kk in xrange(0,n_per):

                
                # add nugget component, apply inverse logit, apply age-correction factor
                chunk = f_chunk_ROW + np.random.normal(loc=0, scale=np.sqrt(V[MCMCrel]), size=f_chunk_ROW.shape)
                chunk = pm.invlogit(chunk.ravel())
                chunk *= facs[np.random.randint(N_facs, size=np.prod(chunk.shape))]
                chunk = chunk.reshape(f_chunk_ROW.shape).squeeze()

                # aggregate through time to obtain spatial-only array for this nugget-realisation
                chunkTMEAN = atleast_2d(np.mean(chunk,-1))
            
                # make a mappng vector for later conversion of arays of this dimension to a vector, and back again
                #ind5km = np.where(chunkTMEAN!=-99999999) 

                # run check that this time-aggregated chunk has same spatial dimensions as time block
                #test=chunk.shape[1:]==chunkTMEAN.shape[0:]
                #if test==False:
                #    print("WARNING !!!!: spatial dimensions of time-aggregated block 'chunkTMEAN' do not match pre-aggregation 'chunk': EXITING!!")
                #    print ('chunk.shape[1:]: '+str(chunk.shape[1:]))
                #    print ('chunkTMEAN.shape[1:]: '+str(chunkTMEAN.shape[1:]))
                #    return(-9999)

                # now expand the 5km PR chunk to match underlying 1km grid (if HiResLowResRatio is 1, will just return copy)
                if(HiResLowResRatio!=1):chunkExpPR = expandGridResPY(chunkTMEAN,HiResLowResRatio)
                if(HiResLowResRatio==1):chunkExpPR = cp.deepcopy(chunkTMEAN)

                # run check that this expanded block has correct dimensions
                test=chunkExpPR.shape==salblim_ROW.shape           
                if test==False:
                    raise RuntimeError ("WARNING !!!!: spatial dimensions of expanded 'chunkExpPR' do not match covariate chunk 'salblim_ROW': EXITING!! ")
                    
                # run check that there are no PR==-9999 pixels (assigned to non-stable pixels in CS code) in stable areas on salblim
                testmatrix = cp.deepcopy(chunkExpPR)
                testmatrix[salblim_ROW == -9999] = 0
                if (np.sum(testmatrix == -9999) > 0):
                    raise RuntimeError ("WARNING!!: ("+str(np.sum(testmatrix== -9999))+") null PR pixels (-9999) found in stable areas in rel "+str(ii)+" , row "+str(jj)+ ": EXITING!!")

                # optionally also derive Ro for this chunk
                if(do_AREALMEANRo|do_POPMEANRo):
                    chunkExpRo=np.zeros(product(shape(chunkExpPR))).reshape(shape(chunkExpPR))
                    nonzero=np.where(chunkExpPR!=0)
                    chunkExpRo[nonzero] = r.PR2Rc(chunkExpPR[nonzero],ES=1, alpha=4.2, c=0.5, s=1, b=0.8, r=355.55/200, gama=0)

                
                #from IPython.Debugger import Pdb
                #Pdb(color_scheme='Linux').set_trace()  

                
                # if needed, also derive PR*pixel area and PR*population version of this chunk and, optionally, equivalent Ro chunks
                if(do_AREALMEANPR):chunkExpPRarea= chunkExpPR * pixarea_ROW 
                if(do_POPMEANPR):chunkExpPRpop= chunkExpPR * grump_ROW
                if(do_AREALMEANRo):chunkExpRoarea= chunkExpRo * pixarea_ROW 
                if(do_POPMEANRo):chunkExpRopop= chunkExpRo * grump_ROW                                 

                if (do_PAR | do_BURDEN):
                
                    if (do_BURDEN):
                        # obtain a burden surface for this chunk as a function of population and PR
                        ## convert PRsurface and POPsurface to vectors before passing, then back=convert afterwards
                        if (HiResLowResRatio!=1): burdenChunk = BurdenPredictorObj.pr5km_pop1km(pr=chunkTMEAN.squeeze(),pop=grump_ROW,pop_pr_res=HiResLowResRatio)
                        if (HiResLowResRatio==1): burdenChunk = BurdenPredictorObj.pr5km_pop5km(pr=chunkTMEAN.squeeze(),pop=grump_ROW)
                
                    # create an ID matrix for this chunk for each endemicity class in each scheme                
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
                            classID[(chunkExpPR>=breaks[cc]) & (chunkExpPR < breaks[cc+1])]=1 
                            #sumClassID=sumClassID+sum(classID)
                            blankarray = {thisbreakname: cp.deepcopy(classID) }

                            # add these blank array to interim PAR dictionary
                            classIDarrays.update(blankarray)

                        # add these sub-dictionaries to PARdict for this scheme
                        classIDarrays = {'classIDarrays':classIDarrays}
                        classIDdict[scheme].update(classIDarrays)

                        # run check that sum of total pixels in all classes in this block match expected total
                        #if sumClassID != np.product(chunkExpPR.shape):
                        #    print "WARNING!! sum of pixels in each class in chunk "+str(jj)+" is "+str(sumClassID)+" != expected total ("+str(np.product(chunkExpPR.shape))+")"

                # loop through each unique country in this chunk: calculate running mean PR,Ro, and running total burden and PAR in each endemicity class in each scheme
                for rr in xrange(0,Nsalb_ROW):
                
                    thiscountry_salbLUT = salbLUT_ROW[rr] 

                    # get this country's ID matrix from countryIDdict dictionary
                    countryID = countryIDdict[str(uniqueSalb_ROW[rr])]

                    # calculate area or population weighted sums of PR and Ro in this country,convert to running mean using known country area or population and add to relevant part of chunk running table
                    if (do_AREALMEANPR):
                        PRsumArea = np.sum(chunkExpPRarea*countryID) 
                        countryAreaMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk] = countryAreaMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk]+(PRsumArea/totalArea[thiscountry_salbLUT])
                        
                    if (do_POPMEANPR):
                        PRsumPop = np.sum(chunkExpPRpop*countryID)
                        countryPopMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk] = countryPopMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk]+(PRsumPop/totalPop[thiscountry_salbLUT])
                        
                    if(do_AREALMEANRo):
                        RosumArea = np.sum(chunkExpRoarea*countryID)
                        countryAreaMeanRorel_ChunkRunning[thiscountry_salbLUT,kk] = countryAreaMeanRorel_ChunkRunning[thiscountry_salbLUT,kk]+(RosumArea/totalArea[thiscountry_salbLUT])
                        
                    if(do_POPMEANRo):
                        RosumPop = np.sum(chunkExpRopop*countryID)
                        countryPopMeanRorel_ChunkRunning[thiscountry_salbLUT,kk] = countryPopMeanRorel_ChunkRunning[thiscountry_salbLUT,kk]+(RosumPop/totalPop[thiscountry_salbLUT])

                    if (do_PAR | do_BURDEN):

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
                                
                                if (do_PAR):
                                    # calculate sum of population in this class and country and add this sum to the relevant part of PARdict_ChunkRunning
                                    PARtemp=grump_ROW * countryClassID
                                    PARsum = np.sum(PARtemp)
                                    PARdict_ChunkRunning[scheme]['PAR'][thisbreakname][thiscountry_salbLUT,kk] = PARdict_ChunkRunning[scheme]['PAR'][thisbreakname][thiscountry_salbLUT,kk] + PARsum

                                    # run check for nonsensical negative PAR or BURDEN value
                                    if PARsum<0.:
                                        print('WARNING!! Negative PAR found - check input population grid. EXITING')
                                        return(-9999)

                                if (do_BURDEN): 
                                    # similarly, calculate sum of burden in this country, and add to relevant part of countryBURDENrel_ChunkRunning 
                                    BURDENtemp = burdenChunk*countryClassID
                                    BURDENsum = np.sum(BURDENtemp)
                                    BURDENdict_ChunkRunning[scheme]['BURDEN'][thisbreakname][thiscountry_salbLUT,kk] = BURDENdict_ChunkRunning[scheme]['BURDEN'][thisbreakname][thiscountry_salbLUT,kk] + BURDENsum
                                        
                                    if BURDENsum<0.:
                                        print('WARNING!! Negative BURDEN found - check input population grid. EXITING')
                                        return(-9999)

        # copy mean PR values for these 500 nugget draws to the main arrays housing all realisations
        if (do_AREALMEANPR): countryAreaMeanPRrel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryAreaMeanPRrel_ChunkRunning
        if (do_POPMEANPR): countryPopMeanPRrel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryPopMeanPRrel_ChunkRunning
        if (do_AREALMEANRo): countryAreaMeanRorel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryAreaMeanRorel_ChunkRunning
        if (do_POPMEANRo): countryPopMeanRorel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryPopMeanRorel_ChunkRunning
        #countryBURDENrel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryBURDENrel_ChunkRunning

        if (do_PAR | do_BURDEN):
            # loop through class schemes
            for ss in xrange(0,Nschemes):
                scheme =  schemeNames[ss]   
                breaknames = breaksDict[scheme]['BREAKNAMES']
                Nclasses=len(breaknames) 

                # .. for each class within each scheme..
                for cc in xrange(0,Nclasses):
                    thisbreakname = breaknames[cc]
                    
                    # copy PAR values for these 500 nugget draws to main arrays housing all realisations
                    if (do_PAR):PARdict[scheme]['PAR'][thisbreakname][:,slice(ii*n_per,(ii*n_per)+n_per,1)] = PARdict_ChunkRunning[scheme]['PAR'][thisbreakname] 

                    # copy BURDEN values for these 500 nugget draws to main arrays housing all realisations
                    if (do_BURDEN):BURDENdict[scheme]['BURDEN'][thisbreakname][:,slice(ii*n_per,(ii*n_per)+n_per,1)] = BURDENdict_ChunkRunning[scheme]['BURDEN'][thisbreakname] 


        
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
    
    # define dictionary to pass to outputDistributedExtractions_country containing some or all of:
    # 1. Array of areal mean PR values per country per realisation (countryAreaMeanPRrel)
    # 2. Array of population mean PR values per country per realisation (countryPopMeanPRrel)
    # 3. Array of areal mean Ro values per country per realisation (countryAreaMeanRorel)
    # 4. Array of population mean Ro values per country per realisation (countryPopMeanRorel)
    # 5. dictionary of burden totals per country per realisation (countryBURDENrel)
    # 6. dictionary of PAR values per country per realisation, along with classification scheme metadata (PARdict)

    returnDict={"startRel":startRelOUT,"endRel":endRelOUT}
    if(do_AREALMEANPR): returnDict.update({"countryAreaMeanPRrel":countryAreaMeanPRrel})
    if(do_POPMEANPR): returnDict.update({"countryPopMeanPRrel":countryPopMeanPRrel})
    if(do_AREALMEANRo): returnDict.update({"countryAreaMeanRorel":countryAreaMeanRorel})
    if(do_POPMEANRo): returnDict.update({"countryPopMeanRorel":countryPopMeanRorel})    
    if(do_PAR): returnDict.update({"PARdict":PARdict})
    if(do_BURDEN): returnDict.update({"BURDENdict":BURDENdict})   
    
    # export extracted arrays 
    outputDistributedExtractions_country(returnDict)
    
    #return(returnDict)
#############################################################################################################################################
def outputDistributedExtractions_country(dict): 

    ''''
    Takes a dictionary which is output from extractSummaries_country and exports 
    .txt files of area/population weighted mean PR or Ro and burden and PAR extractions by country.
    
    Params to pass are:
    
    dict      : output from outputDistributedExtractions_country 
    '''

    # check for error output from extractSummaries_country due to NaNs
    if dict == -9999:
        print "WARNING!! recieved error output from extractSummaries_country() - will not run outputDistributedExtractions_country()" 
        return(-9999)

    # define numpy array of keys in the passed dictionary
    dictKeys = np.array(dict.keys())

    # define which realisations we are dealing with
    startRel=dict['startRel']
    endRel=dict['endRel']

    # construct file suffix indicating realisations in question
    relSuff = '_r'+str(startRel)+'to'+str(endRel)

    # if present, export arrays of area- or population-weighted mean PR or Ro per country per realisation
    if(np.any(dictKeys=='countryAreaMeanPRrel')): np.savetxt(exportPathDistributed_country+'AreaMeanPR_country'+relSuff+'.txt', dict['countryAreaMeanPRrel'])
    if(np.any(dictKeys=='countryAreaMeanRorel')): np.savetxt(exportPathDistributed_country+'AreaMeanRo_country'+relSuff+'.txt', dict['countryAreaMeanRorel'])
    if(np.any(dictKeys=='countryPopMeanPRrel')): np.savetxt(exportPathDistributed_country+'PopMeanPR_country'+relSuff+'.txt', dict['countryPopMeanPRrel'])
    if(np.any(dictKeys=='countryPopMeanRorel')): np.savetxt(exportPathDistributed_country+'PopMeanRo_country'+relSuff+'.txt', dict['countryPopMeanRorel'])    
    
    
    #np.savetxt(exportPathDistributed_country+'BURDEN_country'+relSuff+'.txt', dict['countryBURDENrel'])

    if ((np.any(dictKeys=='PARdict')) | (np.any(dictKeys=='BURDENdict'))):

        Nschemes=len(breaksDict)    
        schemeNames=breaksDict.keys()    

        # loop through classification schemes and clases wihtin them and export array of PAR and/or BURDEN for each
        for ss in xrange(0,Nschemes): 
            scheme=schemeNames[ss]   
            breaknames = breaksDict[scheme]['BREAKNAMES']
            Nclasses=len(breaknames)

            for cc in xrange(0,Nclasses):

                # construct class and scheme suffix
                classSuff = '_'+schemeNames[ss]+'_'+breaknames[cc]

                # export PAR array for this scheme-class
                if(np.any(dictKeys=='PARdict')): np.savetxt(exportPathDistributed_country+'PAR_country'+classSuff+relSuff+'.txt', dict['PARdict'][schemeNames[ss]]['PAR'][breaknames[cc]]) 

                # export BURDEN array for this scheme-class
                if(np.any(dictKeys=='BURDENdict')): np.savetxt(exportPathDistributed_country+'BURDEN_country'+classSuff+relSuff+'.txt', dict['BURDENdict'][schemeNames[ss]]['BURDEN'][breaknames[cc]]) 
#############################################################################################################################################
def extractSummaries_perpixel (slices,a_lo,a_hi,n_per,FileStartRel,FileEndRel,totalN,startRel=None,endRel=None,do_PRMap=False,do_BurdenMap=False,do_RoMap=False):



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


    # define paths to input files according to specified resolution
    if (HiResLowResRatio_PERPIXEL==1):
        salblim_path=salblim5km_path
        salb_path=salb5km_path
        grump_path=grump5km_path
        pixarea_path=pixarea5km_path
        limbnry_path=lim5kmbnry_path
    if (HiResLowResRatio_PERPIXEL==5):
        salblim_path=salblim1km_path
        salb_path=salb1km_path
        grump_path=grump1km_path
        pixarea_path=pixarea1km_path
        limbnry_path=lim1kmbnry_path
    HiResLowResRatio=HiResLowResRatio_PERPIXEL

  
    # check we are actualy asking for one of the map types
    if ((do_PRMap is False) & (do_BurdenMap is False) & (do_RoMap is False)):
        print 'ERROR!!! Not asking for any mapped outputs from extractSummaries_perpixel, EXITING!!'
        return(-9999)


    # construct filepath for this realisation block, and define link
    filename = realisations_path
    filename = filename.replace('FILESTARTREL',str(FileStartRel))
    filename = filename.replace('FILEENDREL',str(FileEndRel))
    hf = tb.openFile(filename)    
    hr = hf.root

    FileStartRel=int(FileStartRel)
    FileEndRel=int(FileEndRel)

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
    N_years = (slices[2].stop - slices[2].start)/12.

    # Get nugget variance and age-correction factors    
    V = hr.PyMCsamples.col('V')[:]    
    facs = mbgw.correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)    

    # if we are extracting burden summaries, import population grid
    if do_BurdenMap==True:
        grump = tb.openFile(grump_path)
        
        # define a function object for later estimation of burden, basedon this grump row (after cnvertig to a vector)
        #ind = np.where(grump5km.root.data[:,:]!=-99999999)
        #POPsurfaceVECTOR=grump5km.root.data[:,:][ind]
        BurdenPredictorObj = BurdenPredictor(hf_name=burdentrace_path, nyr=N_years, burn=0) 

    # define a blank array of zeroes of same size as a single monthly map - that will be duplicated for various uses later
    zeroMap = np.zeros(n_rows*n_cols).reshape(n_rows,n_cols)
    
    # initialise zero matrices that will house running totals
    if do_PRMap==True:
        meanPR = cp.deepcopy(zeroMap)
        meanPR2 = cp.deepcopy(zeroMap)
    if do_BurdenMap==True:
        meanBUR = cp.deepcopy(zeroMap)    
        meanBUR2 = cp.deepcopy(zeroMap)
    if do_RoMap==True:
        meanRo = cp.deepcopy(zeroMap)    
        meanRo2 = cp.deepcopy(zeroMap)

    # initialise dictionary to house probability of class membership running arrays for each scheme/class
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
#        print "realisation "+str(MCMCrel)+" of set "+str(startRel)+" to "+str(endRel)+" ("+str(n_realizations)+" realisations)" 

        #xxx3a = r.Sys_time() 
        # Pull out parasite rate chunk (i.e. import n months of block)    
        tot_slice = (slice(MCMCrel,MCMCrel+1,None),) + slices    
        #f_chunk = hr.realizations[tot_slice][::-1,:,:,:].T   # hr.realizations=[rel,row,col,month]   #f_chunk = [month,col,row,rel]
        #f_chunk = f_chunk[:,:,:,0]                       #f_chunk = [month,col,row]

        # because pyTables is not working properly, manually loop through each month we want to extract and make a ST 3d matrix
        n_months = tot_slice[3].stop - tot_slice[3].start
        f_chunk = np.zeros(1*n_cols*n_rows*n_months).reshape(1,n_rows,n_cols,n_months)
        subsetmonth=0 
        for mm in xrange(tot_slice[3].start,tot_slice[3].stop):
            f_chunk[:,:,:,subsetmonth] = hr.realizations[tot_slice[0],tot_slice[1],tot_slice[2],mm]
            subsetmonth=subsetmonth+1
        #f_chunk = f_chunk[::-1,:,::-1,:].T[:,:,:,0]
        f_chunk = f_chunk.squeeze()

        ########TEMP###########
        #set missing vlaues in f block to 0
        #print(sum(np.isnan(f_chunk)))
        #f_chunk[np.isnan(f_chunk)]=0
        #print(sum(np.isnan(f_chunk)))
        ####################################

        print('mean of f_chunk: '+str(np.mean(f_chunk)))
         
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
            chunk = f_chunk + np.random.normal(loc=0, scale=np.sqrt(V[MCMCrel]), size=f_chunk.shape)
            chunk = pm.invlogit(chunk.ravel())
            chunk *= facs[np.random.randint(N_facs, size=np.prod(chunk.shape))]
            chunk = chunk.reshape(f_chunk.shape).squeeze()

#            # if we are working out Ro, need to calculate before temporal aggregation of PR cube
#            if do_RoMap==True:
#                chunk_Ro=np.zeros(product(shape(chunk))).reshape(shape(chunk))
#                nonzero=np.where(chunk!=0)
#                chunk_Ro[nonzero] = r.PR2Rc(chunk[nonzero],ES=1, alpha=4.2, c=0.5, s=1, b=0.8, r=355.55/200, gama=0)

            # aggregate through time to obtain spatial-only array for this nugget-realisation
            chunkTMEAN = np.atleast_2d(np.mean(chunk,-1))
            print('mean of chunkTMEAN: '+str(np.mean(chunkTMEAN)))

            if do_RoMap==True:
#                chunkTMEAN_Ro = np.atleast_2d(np.mean(chunk_Ro,-1))
#                print('mean of chunkTMEAN_Ro: '+str(np.mean(chunkTMEAN_Ro)))
                chunkTMEAN_Ro=np.zeros(product(shape(chunkTMEAN))).reshape(shape(chunkTMEAN))
                nonzero=np.where(chunkTMEAN!=0)
                chunkTMEAN_Ro[nonzero] = r.PR2Rc(chunkTMEAN[nonzero],ES=1, alpha=4.2, c=0.5, s=1, b=0.8, r=355.55/200, gama=0)

#            from IPython.Debugger import Pdb
#            Pdb(color_scheme='Linux').set_trace()                   
            
            # increment running mean matrices 
            if do_PRMap==True:
                meanPR = meanPR + (chunkTMEAN/totalN)
                meanPR2 = meanPR2 + (np.square(chunkTMEAN)/totalN)

            if do_RoMap==True:
                meanRo = meanRo + (chunkTMEAN_Ro/totalN)
                meanRo2 = meanRo2 + (np.square(chunkTMEAN_Ro)/totalN)

            if do_BurdenMap==True:
                burdenChunk = BurdenPredictorObj.pr5km_pop5km(pr=chunkTMEAN,pop=grump.root.data[:,:])
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

    print('mean of meanPR before return as .gz: '+str(np.mean(meanPR)))

    ## export running meanPR and meanPR2 array
    if do_PRMap==True:
        np.savetxt(exportPathDistributed_perpixel+"meanPR_perpixel"+relSuff+".gz",meanPR)
        np.savetxt(exportPathDistributed_perpixel+"meanPR2_perpixel"+relSuff+".gz",meanPR2)

    ## export running meanRo and meanRo2 array
    if do_RoMap==True:
        np.savetxt(exportPathDistributed_perpixel+"meanRo_perpixel"+relSuff+".gz",meanRo)
        np.savetxt(exportPathDistributed_perpixel+"meanRo2_perpixel"+relSuff+".gz",meanRo2)

    ## export running meanBUR and meanBUR2 array
    if do_BurdenMap==True:
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
#############################################################################################################################################def queryRealizationsInFolder(folderPath,relPath,VERBOSE=True):

    '''
    Queries a specified bucket for realisation files and returns dictionary with N realisations, and lists of start and end realisation numbers

    params to pass:
    folderPath     : (string) path to folder containing realisations
    relPath       : (string) a generic filename for the realisatios in this bucket (same format as that in extract params, with FILESTARTREL
                    and FILEENDREL suffixes denoting the firs and last realisatoin couintained in each realisation filepath to target directory. If this includes new directories, these will be built if possible
    '''

    # define reference filepath for later use - removing the .hdf5 suffix                   
    TEMPrelPath = relPath.rpartition('.')[-3]

    # check folder exists
    if (checkAndBuildPaths (folderPath,VERBOSE=True,BUILD=False) ==-9999):
        print 'WARNING!!! requested folder "'+str(folderPath)+'" does not exist !!!'
        return(-9999)

    # query folder and define list of filenames therein
    filelist =  os.listdir(folderPath)
    NfilesInFolder = len(filelist)

    Nrealisations = 0
    Nfiles = 0
    StartRelList = []
    EndRelList = []

    # loop through filenames in bucket, check they are a realisatoin file, and if so extract start and end realistaion info
    for i in xrange(0,NfilesInFolder):

        # get the ith filename in the bucket
        fileN=filelist[i]
        
        # first check it is an hdf5 file
        if fileN.rpartition('.')[-1]=='hdf5':
            
            # truncate this filename to remove extension
            TEMPPath = fileN.rpartition('.')[-3]
        
            # compare this filename to the generic realisation filename (after removing the realisation-specific suffixes), if they match then assume this is a realisation file
            if (TEMPPath.split('_')[0:-2:1]==TEMPrelPath.split('_')[0:-2:1]):
                
                # extract start and end realisation numbers from file suffixes
                StartRelList.append(TEMPPath.split('_')[-2])
                EndRelList.append(TEMPPath.split('_')[-1])
                
                # increment counter of valid realisation files in this bucket
                Nfiles+=1
                
                # and increment counter of nuber of individual realisations in this bucket (a file may have more than one realisation)
                Nrealisations=Nrealisations +  int(TEMPPath.split('_')[-1]) - int(TEMPPath.split('_')[-2])                     

    # return dictionary with start and end realisation lists and number of realisations
    returnDict = {"Nfiles":Nfiles,"Nrealisations":Nrealisations,"StartRelList":StartRelList,"EndRelList":EndRelList}
    return returnDict
#############################################################################################################################################
