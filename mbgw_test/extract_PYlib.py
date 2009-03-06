# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################


# import python libraries
from rpy import *
import numpy as np
import copy as cp
import pymc as pm
import mbgw
import tables as tb
import pylab as pl
import matplotlib
matplotlib.interactive(True)
from mbgw.joint_simulation import *
from mbgw import st_cov_fun

# import R functions
r.source('extract_Rlib.R')
plotmonthPY = r['plotmonth']
getTimeMeanPY = r['getTimeMean']
expandGridResPY=r['expandGridRes']

# import parameters from param file
from extract_params import *

# check filepaths stated in parameter file
checkAndBuildPaths(filename,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPath,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathCombined,VERBOSE=True,BUILD=True)
checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=True)


#############################################################################################################################################
def examineSalb (salblim1km_path,uniqueSalb_path,pixelN_path):
    salblim1km = tb.openFile(salblim1km_path, mode = "r")
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

    uniqueSalb.tofile(uniqueSalb_path,sep=",")
    count.tofile(pixelN_path,sep=",")

#############################################################################################################################################
def extractSTaggregations (slices,a_lo,a_hi,n_per,startRel,endRel):

    ####TEMP
    #XXXa=r.Sys_time()
    #xxx1=xxx2=xxx3=xxx4=xxx5=xxx6=xxx7=xxx8=xxx9=xxx10=xxx11=xxx12=xxx13=xxx14=xxx15=xxx16=xxx17=xxx18=xxx19=xxx20=0
    ########
    
    #xxx1a = r.Sys_time()
    # define basic parameters
    slices = tuple(slices)     
    hf = tb.openFile(filename)    
    hr = hf.root    
    #n_realizations = len(hr.realizations)    
    n_realizations = (endRel - startRel)+1
    n_rows=len(hr.lat_axis)
    N_facs = int(1e5)    

    # Get nugget variance and age-correction factors    
    V = hr.PyMCsamples.col('V')[:]    
    facs = mbgw.correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)    

    # open link to salb grid (masked to stable areas only) and population grid    
    salblim1km = tb.openFile(salblim1km_path, mode = "r")    
    gr001km = tb.openFile(gr001km_path, mode = "r")    

    # perform check that the number of rows and columns is the same in both 1km grids
    if len(salblim1km.root.lat) != len(gr001km.root.lat):
        print 'WARNING!! 1km row numbers do not correspond: salblim1km has '+str(len(salblim1km.root.lat))+' and gr001km has '+str(len(gr001km.root.lat))
    if len(salblim1km.root.long) != HiResLowResRatio*len(hr.lon_axis):
        print 'WARNING!! col numbers do not correspond: salblim1km has '+str(len(salblim1km.root.long))+' and gr001km has '+str(len(gr001km.root.long))

    # perform check that the number of rows and columns is in the correct ratio to those of input 5km grid
    if len(salblim1km.root.lat) != HiResLowResRatio*len(hr.lat_axis):
        print 'WARNING!! 1km and 5km row numbers do not correspond: salblim1km has '+str(len(salblim1km.root.lat))+' and 5km rows * HiResLowResRatio is '+str(HiResLowResRatio*len(hr.lat_axis))
    if len(salblim1km.root.long) != HiResLowResRatio*len(hr.lon_axis):
        print 'WARNING!! 1km and 5km col numbers do not correspond: salblim1km has '+str(len(salblim1km.root.long))+' and 5km cols * HiResLowResRatio is '+str(HiResLowResRatio*len(hr.long_axis))

    # get list of unique salb IDs and count of pixels in each..
    # ..first check that Salb grid has been pre-examined using examineSalb and lists of unique IDs and N pixels exist, if not then re-run examineSalb
    try:
        uniqueSalb=fromfile(uniqueSalb_path,sep=",")
        pixelN=fromfile(pixelN_path,sep=",")
    except IOError:
        print 'WARNING!! files '+pixelN_path+" or "+uniqueSalb_path+" not found: running examineSalb"
        examineSalb (salblim1km_path,uniqueSalb_path,pixelN_path)

    uniqueSalb=fromfile(uniqueSalb_path,sep=",")    
    pixelN=fromfile(pixelN_path,sep=",")        
    Nsalb=len(uniqueSalb)    

    # intialise empty arrays (e.g. 87 countries * N realisations) for mean PR and PAR in each class of each scheme..    
    countryMeanPRrel = repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per)     

    # intialise empty arrays (e.g. 87 countries * N realisations) for PAR in each class of each scheme..housed in PAR dictionary PARdict
    Nschemes=len(breaksDict)    
    schemeNames=breaksDict.keys()    
    PARdict=cp.deepcopy(breaksDict)
    #xxx1 = xxx1 + (r.Sys_time() - xxx1a)
    
    #xxx2a = r.Sys_time()
    # ..loop through each classification scheme 
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = PARdict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames) 

        # define additional sub-dictionary to add to PARdict to house arrays for PAR per class per scheme per country
        PAR = {}

        # .. for each class within each scheme..
        for cc in xrange(0,Nclasses):
            thisbreakname = breaknames[cc]
            
            # define an empty array for this scheme-class to house PAR per country realisations for each class
            blankarray = {thisbreakname: repeat(None,n_realizations*n_per*Nsalb).reshape(Nsalb,n_realizations*n_per) }

            # add this blank array to interim PAR dictionary
            PAR.update(blankarray)

        # add this sub-dictionary to PARdict for this scheme
        PAR = {'PAR':PAR}
        PARdict[scheme].update(PAR)
    #xxx2 = xxx2 + (r.Sys_time() - xxx2a)

    # loop through each realisation
    for ii in xrange(0,n_realizations): #1:500 realisations n_realizations   
    
        # define which realisatin this relates to in global set from MCMC
        MCMCrel = startRel+ii 

        print "realisation "+str(MCMCrel)+" of set "+str(startRel)+" to "+str(endRel)+" ("+str(n_realizations)+" realisations)"    

        #xxx3a = r.Sys_time() 
        # Pull out parasite rate chunk (i.e. import n months of block)    
        tot_slice = (slice(MCMCrel,MCMCrel+1,None),) + slices    
        f_chunk = hr.realizations[tot_slice][::-1,:].T   # hr.realizations=[rel,row,col,month]   #f_chunk = [month,col,row,rel]
        f_chunk = f_chunk[:,:,:,0]                       #f_chunk = [month,col,row]

        ########TEMP###########
        #set missing vlaues in f block to 0
        #from scipy.io import write_array
        #write_array('/home/pwg/MBGWorld/extraction/temp_PRrel1.txt', f_chunk[0,:,:])
        #print(sum(isnan(f_chunk)))
        f_chunk[isnan(f_chunk)]=0
        #print(sum(isnan(f_chunk)))
        ####################################
        #xxx3 = xxx3 + (r.Sys_time() - xxx3a)

        #xxx4a = r.Sys_time() 
        # run check that there are no missing values in this f chunk
        if sum(isnan(f_chunk))>0:
            print "WARNING!! found "+str(sum(isnan(f_chunk)))+" NaN's in realisation "+str(MCMCrel)+" EXITING!!!"
            return()

        # initialise arrays to house running mean PR whilst we loop through chunks and nugget draws..
        countryMeanPRrel_ChunkRunning = repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per)

        # initialise corresponding arrays for running PAR (housed in PR dict, so simply ensure reset to zero for this realisation)..
        PARdict_ChunkRunning=cp.deepcopy(breaksDict)

        # ..loop through each classification scheme.. 
        for ss in xrange(0,Nschemes):
            scheme = schemeNames[ss]    
            breaknames = PARdict_ChunkRunning[scheme]['BREAKNAMES']
            Nclasses=len(breaknames) 

            # define additional sub-dictionaries to add to PARdict_ChunkRunning to house temporary arrays
            PAR = {}

            # .. for each  class within each scheme..
            for cc in xrange(0,Nclasses):
                thisbreakname = breaknames[cc]
            
                # define an empty array for this scheme-class to house running PAR sum
                blankarray = {thisbreakname: repeat(0.,n_per*Nsalb).reshape(Nsalb,n_per) }

                # add these blank array to interim PAR dictionary
                PAR.update(blankarray)

            # add these sub-dictionaries to PARdict for this scheme
            PAR = {'PAR':PAR}
            PARdict_ChunkRunning[scheme].update(PAR)
        #xxx4 = xxx4 + (r.Sys_time() - xxx4a)
        
        # loop through each row (or multiple rows in a slice) of 5km realisation grid..
        #n_slices = n_rows/rowsInslice5km

        interimCnt=0 
        for jj in xrange(0,n_rows): 

            interimCnt=interimCnt+1
            if interimCnt==10:
                print('    on slice '+str(jj)+' of '+str(n_rows))
                interimCnt=0
                    
            #xxx5a = r.Sys_time() 
            # get row of 5km PR surface accross all months in chunk  (assumes f_chunk is correct way up i.e. map view)
            f_chunk_ROW = f_chunk[:,jj,:]

            # get corresponding 5 rows of 1km Salb and population surface (assumes they are correct way up i.e. map view)
            startRow1km=jj*HiResLowResRatio
            endRow1km=startRow1km+(HiResLowResRatio)
            salblim1km_ROW = salblim1km.root.data[slice(startRow1km,endRow1km,1),:]
            gr001km_ROW = gr001km.root.data[slice(startRow1km,endRow1km,1),:]
            
            # define a blank array of zeroes of same size as 1km chunk - that will be duplicated for various uses later
            zeroChunk = zeros(product(gr001km_ROW.shape)).reshape(gr001km_ROW.shape)
            #xxx5 = xxx5 + (r.Sys_time() - xxx5a)

            #plotMapPY(salblim1km.root.data[:,:],NODATA=-9999)
            #plotMapPY(salblim1km.root.data[slice(0,100,1),:],NODATA=-9999)
            #plotMapPY(f_chunk[0,:,:])
            #plotMapPY(f_chunk[0,slice(0,100,1),:])

            #xxx6a = r.Sys_time() 
            # how many unique salb IDs in these rows (after removing -9999 cells from e.g. sea)?
            uniqueSalb_ROW = unique(salblim1km_ROW)
            #uniqueSalb_ROW = uniqueSalb_ROW[uniqueSalb_ROW!=-9999]
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
                sumCountryID = sumCountryID + sum(countryIDmatrix)
                tmpdict={str(uniqueSalb_ROW[rr]):cp.deepcopy(countryIDmatrix) }
                countryIDdict.update(tmpdict)
            #xxx7 = xxx7 + (r.Sys_time() - xxx7a)
            
            # run check that sum of total pixels in all countries in this block match expected total
            if sumCountryID != product(salblim1km_ROW.shape):
                print "WARNING!! sum of pixels in each country in chunk "+str(jj)+" is "+str(sumCountryID)+" != expected total ("+str(product(salblim1km_ROW.shape))+")"

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
                #chunkTMEAN = array([getTimeMeanPY(chunk,TimeDim=1)])
                #print(shape(chunkTMEAN))
                #print(mean(chunkTMEAN))
                #print(type(chunkTMEAN))
                chunkTMEAN = atleast_2d(np.mean(chunk,0))
                #print(shape(chunkTMEAN))
                #print(mean(chunkTMEAN))
                #print(type(chunkTMEAN))
                #print(" ")
                #xxx9 = xxx9 + (r.Sys_time() - xxx9a)

                # run check that this time-aggregated chunk has same spatial dimensions as time block
                test=chunk.shape[1:]==chunkTMEAN.shape[1:]
                if test==False:
                    print("WARNING !!!!: spatial dimensions of time-aggregated block 'chunkTMEAN' do not match pre-aggregation 'chunk'")

                #xxx10a = r.Sys_time()
                # now expand the 5km PR chunk to match underlying 1km grid
                chunkExp = expandGridResPY(chunkTMEAN,HiResLowResRatio)
                #xxx10 = xxx10 + (r.Sys_time() - xxx10a)                

                # run check that this expanded block has correct dimensions
                test=chunkExp.shape==salblim1km_ROW.shape           
                if test==False:
                    print("WARNING !!!!: spatial dimensions of expanded 5km 'chunkExp' do not match 1km covariate chunk 'salblim1km_ROW' ")            

                # create an ID matrix for this chunk for each class in each scheme                
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
                    #if sumClassID != product(chunkExp.shape):
                    #    print "WARNING!! sum of pixels in each class in chunk "+str(jj)+" is "+str(sumClassID)+" != expected total ("+str(product(chunkExp.shape))+")"
                #xxx11 = xxx11 + (r.Sys_time() - xxx11a)

                # loop through each unique country in this chunk: calculate running mean PR and PAR
                for rr in xrange(0,Nsalb_ROW):
                
                    thiscountry_salbLUT = salbLUT_ROW[rr] 

                    # get this country's ID matrix from countryIDdict dictionary
                    countryID = countryIDdict[str(uniqueSalb_ROW[rr])]

                    #xxx12a = r.Sys_time()
                    # calculate sum of PR in this country,convert to running mean  using known conutry pixel count, and add to the relevant part of countryMeanPRrel_ChunkRunning
                    PRsum = sum(chunkExp*countryID)
                    countryMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk] = countryMeanPRrel_ChunkRunning[thiscountry_salbLUT,kk]+(PRsum/pixelN[thiscountry_salbLUT])
                    #xxx12 = xxx12 + (r.Sys_time() - xxx12a)

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
                            #
                            countryClassID = countryID*classID   
                            #xxx13 = xxx13 + (r.Sys_time() - xxx13a)
                            
                            #xxx14a = r.Sys_time()
                            # calculate sum of population in this class and country and add this sum to the relevant part of PARdict_ChunkRunning
                            PARtemp=gr001km_ROW * countryClassID
                            #xxx14 = xxx14 + (r.Sys_time() - xxx14a)
                            #xxx15a = r.Sys_time()
                            PARsum = sum(PARtemp)
                            #xxx15 = xxx15 + (r.Sys_time() - xxx15a)
                            
                            #xxx16a = r.Sys_time()
                            PARdict_ChunkRunning[scheme]['PAR'][thisbreakname][thiscountry_salbLUT,kk] = PARdict_ChunkRunning[scheme]['PAR'][thisbreakname][thiscountry_salbLUT,kk] + PARsum
                            #xxx16 = xxx16 + (r.Sys_time() - xxx16a)

                            # run check for nonsensical negative PAR value
                            if PARsum<0.:
                                print('WARNING!! Negative PAR found - check input population grid. EXITING')
                                return()
                    #xxx13 = xxx13 + (r.Sys_time() - xxx13a)

        # copy mean PR values for these 500 nugget draws to the main arrays housing all realisations
        #xxx17a = r.Sys_time()
        countryMeanPRrel[:,slice(ii*n_per,(ii*n_per)+n_per,1)] = countryMeanPRrel_ChunkRunning

        # loop through class schemes
        for ss in xrange(0,Nschemes):
            scheme =  schemeNames[ss]   
            breaknames = PARdict[scheme]['BREAKNAMES']
            Nclasses=len(breaknames) 

            # .. for each class within each scheme..
            for cc in xrange(0,Nclasses):
                thisbreakname = breaknames[cc]
                
                # copy PAR values for these 500 nugget draws to main arrays housing all realisations
                PARdict[scheme]['PAR'][thisbreakname][:,slice(ii*n_per,(ii*n_per)+n_per,1)] = PARdict_ChunkRunning[scheme]['PAR'][thisbreakname] 
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

    # return dictionary containing:
    # 1. Array of mean PR values per country per realisation (countryMeanPRrel)
    # 2. dictionary of PAR values per country per realisation, along with classification scheme metadata (PARdict)

    returnDict={"startRel":startRel,"endRel":endRel,"countryMeanPRrel":countryMeanPRrel,"PARdict":PARdict}
    return(returnDict)
 
#############################################################################################################################################
def outputExtraction(dict): 
    
    from scipy.io import write_array
    
    # define which realisations we are dealing with
    startRel=dict['startRel']
    endRel=dict['endRel']
    
    # construct file suffix indicating realisations in question
    relSuff = '_r'+str(startRel)+'to'+str(endRel)
    
    # export array of mean PR per country per realisation
    write_array(exportPath+'meanPR'+relSuff+'.txt', dict['countryMeanPRrel'])
    
    # loop through classification schemes and clases wihtin them and export array of PAR for each
    schemes = dict['PARdict'].keys()    
    Nschemes = len(schemes)
    
    for ss in xrange(0,Nschemes):
        classes = dict['PARdict'][schemes[ss]]['PAR'].keys()
        Nclasses = len(classes)
        
        for cc in xrange(0,Nclasses):
            
            # construct class and scheme suffix
            classSuff = '_'+schemes[ss]+'_'+classes[cc]

            # export PAR array for this scheme-class
            write_array(exportPath+'PAR'+classSuff+relSuff+'.txt', dict['PARdict'][schemes[ss]]['PAR'][classes[cc]]) 
#############################################################################################################################################

a=r.Sys_time()
ExtractedDict = extractSTaggregations([slice(None,None,None), slice(None,None,None), slice(0,12,None)],2,10,1,0,1)
print("TOTAL TIME: "+(str(r.Sys_time()-a)))
outputExtraction(ExtractedDict) 