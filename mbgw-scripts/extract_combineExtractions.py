# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################

# import python libraries
import os
import numpy as np
import copy as cp
import tables as tb
from extract_PYlib import *
from rpy import *
from map_utils import quantile_funs as qs
from map_utils import getAsciiheaderFromTemplateHDF5
from map_utils import exportAscii
from UserDict import UserDict

# import parameters
from extract_params import *

# import some r functions
r.source(utilFolder+'GeneralUtility.R')
writeTableWithNamesPY = r['writeTableWithNames'] 

#############################################################################################################################################
def FnamesByVariableAndOrderedByRel(fnames,variableName=None):

    FnamesVariable = []
    FnamesVariable_StartRel = []

    for fname in fnames:
        nameparts = deconstructFilename(fname)
        if variableName is None:
            FnamesVariable.append(fname)
            FnamesVariable_StartRel.append(nameparts['startRel'])        
        
        if variableName is not None:
            if nameparts['variable']==variableName:
                FnamesVariable.append(fname)
                FnamesVariable_StartRel.append(nameparts['startRel'])

    FnamesVariableOrdered = []

    while len(FnamesVariable_StartRel) >0:
        
        # determine position of smallest startRel
        minelement = np.where(FnamesVariable_StartRel==np.min(FnamesVariable_StartRel))
        minelement=minelement[0][0]
        minelement=int(minelement)

        # pop out this element from both lists
        null = FnamesVariable_StartRel.pop(minelement)
        FnamesVariableOrdered.append(FnamesVariable.pop(minelement))

    return(FnamesVariableOrdered)
#############################################################################################################################################
def deconstructFilename (fname):
    # deconstruct filename to assign variables relating to this file 
    part1 = fname.partition('_')
    variable = part1[0]

    if ((variable==str('PAR')) | (variable==str('PCM'))):
       part2 = part1[2].partition('_')
       extractiontype = part2[0]
       part3 = part2[2].partition('_')
       scheme = part3[0]
       part4 = part3[2].partition('_')
       breakname = part4[0]
       part5=part4[2].partition('r')
       part6=part5[2].partition('to')
       startRel = int(part6[0])
       part7=part6[2].partition('.')
       endRel = int(part7[0])
       returnDict={'variable':variable,'extractiontype':extractiontype,'scheme':scheme,'breakname':breakname,'startRel':startRel,'endRel':endRel}
       return(returnDict)

    if ((variable==str('meanPR')) | (variable==str('BURDEN')) | (variable==str('meanPR2')) | (variable==str('meanBUR')) | (variable==str('meanBUR2'))):       
       part2= part1[2].partition('_')
       extractiontype = part2[0]        
       part3=part2[2].partition('r')
       part4=part3[2].partition('to')
       startRel = int(part4[0])
       part5=part4[2].partition('.')
       endRel = int(part5[0])
       returnList=variable,extractiontype,startRel,endRel
       returnDict={'variable':variable,'extractiontype':extractiontype,'startRel':startRel,'endRel':endRel}
       return(returnDict)
       
    print "WARNING!!! file "+fname+" does not contain expected variables"
    return(-9999)

#############################################################################################################################################
def copySubTableToMain(subfname,maintable,n_per,Nsalb,filledCols):

    # read in file
    inputTable = np.loadtxt(exportPathDistributed_country+subfname)

    # if there is only a single column of input (which would mean 1 nugget realisation from 1 full realisation), then deal with dimensionality
    if len(shape(inputTable)) == 1:
        inputTable = np.atleast_2d(inputTable).T

    ## query attributes of input table
    name_parts = deconstructFilename(subfname)    
    startRel = name_parts['startRel']    
    endRel = name_parts['endRel']    
    ncols = inputTable.shape[1]    
    nrel = (endRel-startRel)    

    # check attributes are as expected    
    if inputTable.shape[0] != Nsalb:    
        print( "Warning!!! nrows of "+subfname+" is "+str(inputTable.shape[0])+" != expected "+str(Nsalb))    
    if ncols/nrel != n_per:    
        print( "Warning!!! ncols of "+subfname+" is "+str(ncols/nrel)+" != expected "+str(n_per))    

    # copy contents of this file to the global array in the relevant position    
    #startCol = startRel*n_per    
    #endCol = startCol + (nrel*n_per)
    startCol = filledCols    
    endCol = filledCols + (nrel*n_per)

    print str(subfname)
    print 'startCol: '+str(startCol)
    print 'endCol: '+str(endCol)
    print 'filledCols: '+str(filledCols)
    print '\n'

    maintable[:,startCol:endCol:1]=inputTable
    
    # increment number of columns now filled on maintable
    filledCols = endCol

    # return modified main table and number of columns already filled
    returnDict = {'maintable':maintable,'filledCols':filledCols}
    return returnDict
#############################################################################################################################################
def makeGlobalArray_contVariables(variableName,blankarray,n_per,Nsalb):

    globalarray = cp.deepcopy(blankarray)
    filledCols = 0

    # get sub list of filename in directory ordered by start ralisation, and optionally containng the variable name
    allfnames = os.listdir(exportPathDistributed_country)
    fnames = FnamesByVariableAndOrderedByRel(allfnames,variableName)

    for fname in fnames:
        if fname.find(variableName+'_')!=-1:

            # copy this array to correct position on global array
            returnDict =  copySubTableToMain(fname,globalarray,n_per,Nsalb,filledCols)
            globalarray = returnDict['maintable']
            filledCols = returnDict['filledCols']

    np.savetxt(exportPathCombined_country+variableName+".txt",globalarray)

    # optionally also create a summary table - mean,SD quantiles etc over all realisations per country
    if summaryStats!=None:
        summObj = getSummariesPerCountry(globalarray)
        summTable = summObj[0]
        summNames = summObj[1]
        writeTableWithNamesPY(summTable,names=summNames,filepath=exportPathCombined_country+variableName+"_summary.csv",SEP=",")

    return(globalarray) 

#############################################################################################################################################
def makeGlobalArray_categoVariables(variableName,blankarray,n_per,Nsalb):

    Nschemes=len(breaksDict)    
    schemeNames=breaksDict.keys()    

    # ..loop through each classification scheme 
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = breaksDict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames)

        # .. for each class within each scheme..
        for cc in xrange(0,Nclasses):
            thisbreakname = breaknames[cc]

            # loop through all meanPR_* files and add them to global array, then export global array
            globalarray = cp.deepcopy(blankarray)

            # get sub list of filename in directory ordered by start ralisation, and optionally containng the variable name
            allfnames = os.listdir(exportPathDistributed_country)
            fnames = FnamesByVariableAndOrderedByRel(allfnames,variableName)
            print fnames
            filledCols = 0
            for fname in fnames:
                print fname
                if (fname.find(variableName+'_')!=-1) & (fname.find(scheme)!=-1) & (fname.find(thisbreakname)!=-1) :

                    # copy this array to correct position on global array
                    returnDict =  copySubTableToMain(fname,globalarray,n_per,Nsalb,filledCols)
                    globalarray = returnDict['maintable']
                    filledCols = returnDict['filledCols']

            np.savetxt(exportPathCombined_country+variableName+'_'+scheme+'_'+thisbreakname+'.txt',globalarray)

            # optionally also create a summary table - mean,SD quantiles etc over all realisations per country
            if summaryStats!=None:
                summObj = getSummariesPerCountry(globalarray)
                summTable = summObj[0]
                summNames = summObj[1]
                writeTableWithNamesPY(summTable,names=summNames,filepath=exportPathCombined_country+variableName+'_'+scheme+'_'+thisbreakname+"_summary.csv",SEP=",")

#############################################################################################################################################
def getSummariesPerCountry(globalarray):

    # ..first check that Salb grid has been pre-examined using examineSalb and lists of unique IDs and N pixels exist, if not then re-run examineSalb
    try:
        uniqueSalb=fromfile(uniqueSalb_path,sep=",")
        pixelN=fromfile(pixelN_path,sep=",")
    except IOError:
        print 'WARNING!! files '+pixelN_path+" or "+uniqueSalb_path+" not found: running examineSalb"
        temp=examineSalb (salblim1km_path,ignore=np.array([-9999]))
        uniqueSalb=temp['uniqueSalb']
        pixelN=temp['count'] 

    outTable=atleast_2d(uniqueSalb).T
    colNames = list(['salbID'])

    #print outTable

    statArray=(None)
    for stat in summaryStats:
        #print stat
        if stat=='mean':
            outTable = np.append(outTable,atleast_2d(np.mean(globalarray,1)).T,axis=1)
            colNames.extend(['mean'])
            #print outTable
        if stat=='SD':
            #print atleast_2d(np.std(globalarray,1)).T 
            outTable = np.append(outTable,atleast_2d(np.std(globalarray,1)).T,axis=1)
            colNames.extend(['SD'])
            #print outTable
        if stat=="quantiles":
            #print qs.row_quantile(globalarray, summaryStats[stat])
            outTable = np.append(outTable,qs.row_quantile(globalarray, summaryStats[stat]),axis=1)
            colNames.extend(r.as_character(summaryStats[stat]))        
            #print outTable

#    print outTable
#    print colNames 
    return outTable,colNames

#############################################################################################################################################
def combineDistribExtractions_country():

    # loop through all files in 'exportPathDistributed_country' directory, get list of unique realization ranges, and thereby
    # calculate total top-level realizations present accross output files in directory
    relRanges=list()  
    i=0
    for fname in os.listdir(exportPathDistributed_country):

        # does this file contain the '_r' string and is it a text file - check that this is the extraction output (although should not be anything else in directory)
        if fname.find('_r') == -1 | fname.find('.txt') == -1: continue

        # record name of first file for later use
        if i==0:
            firstfname = fname
            i=1

        # deconstruct filename
        name_parts = deconstructFilename (fname)
        relRange=name_parts['startRel'],name_parts['endRel']
        relRanges.append(relRange)

    ## get array of unique rel ranges
    uniqueRelRanges=np.unique(relRanges)
    Nunique = uniqueRelRanges.shape[0] 

    ## calculate total number of top-level realisations present accross these files
    n_realizations_infiles=0
    for i in xrange(0,Nunique):
        subtotal = (uniqueRelRanges[i][1] - uniqueRelRanges[i][0])
        n_realizations_infiles = n_realizations_infiles + subtotal

    # import first file and obtain dimensions to reveal N_uniqueSalb and n_per 
    name_parts = deconstructFilename (firstfname)
    firstfile=np.loadtxt(exportPathDistributed_country+firstfname)
 
    # if there is only a single column of input (which would mean 1 nugget realisation from 1 full realisation), then deal with dimensionality
    if len(shape(firstfile)) == 1:
        firstfile = np.atleast_2d(firstfile).T

    tableshape = firstfile.shape
    Nsalb = tableshape[0]
    ncols = tableshape[1] 
    n_per = ncols/(name_parts['endRel']-name_parts['startRel'])

    # define generic blank array to house combined tables for each variable
    blankarray = np.repeat(-9999.,n_realizations_infiles*n_per*Nsalb).reshape(Nsalb,n_per*n_realizations_infiles)

    # loop through all meanPR_* files and add them to global array, then export global array, and optionallly accopmanying summary table
    temp=makeGlobalArray_contVariables('meanPR',blankarray,n_per,Nsalb)

    # loop through all BURDEN_* files and add them to global array, then export global array, and optionallly accopmanying summary table
    temp=makeGlobalArray_contVariables('BURDEN',blankarray,n_per,Nsalb)

    # loop through all PAR_* files and add them to global arrays, then export global arrays, and optionallly accopmanying summary table
    makeGlobalArray_categoVariables('PAR',blankarray,n_per,Nsalb)

    return()
#############################################################################################################################################
def combineDistribExtractions_perpixel():

    # import limits mask to neaten up and provide ascii template for final output arrays
    mask = tb.openFile(lim5kmbnry_path)
    hdrDict = getAsciiheaderFromTemplateHDF5(lim5kmbnry_path)

    # import first file to establish array dimensions (should all be the same in this directory, and all subsequent imports will be checked for consistency)
    dirList = os.listdir(exportPathDistributed_perpixel)
    for i in  xrange(0,len(dirList)):
        fname = dirList[i]
        if fname.find('_r') == -1 | fname.find('.gz') == -1:
            print 'WARNING!! file '+str(fname)+' in '+str(exportPathDistributed_perpixel)+') is not of correct type, trying next file!!!'
            if (i==(len(dirList)-1)):
                print 'ERROR!!! no suitable files found in '+str(exportPathDistributed_perpixel)+' : EXITING!!!'
                return(-9999)
            continue
        else:
            # import first suitable file, compare its shape to mask grid, and if OK use as reference shape
            referenceshape = np.loadtxt(exportPathDistributed_perpixel+fname).shape
            if referenceshape!=mask.root.data[:,:].shape:
                print 'ERROR!!! mask shape is '+str(mask.root.data[:,:].shape)+ 'but file '+str(fname)+' has shape '+str(referenceshape)+' : EXITING!!!!'
                return(-9999)
            break

    # establish zero array in this shape that will be duplicated
    zeroMap = zeros(product(referenceshape)).reshape(referenceshape)
     
    # initialise zero arrays to sum over for means, and counters for checking 
    meanPR = cp.deepcopy(zeroMap)
    meanPR2 = cp.deepcopy(zeroMap) 
    meanBUR = cp.deepcopy(zeroMap) 
    meanBUR2 = cp.deepcopy(zeroMap)

    meanPRtally = meanPR2tally = meanBURtally = meanBUR2tally = 0 

    # initialise dictionary of PCM arrays to sum over for each scheme/class
    Nschemes=len(breaksDict)    
    schemeNames=breaksDict.keys()    
    PCMdict=cp.deepcopy(breaksDict)
    
    ## ..loop through each classification scheme 
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = PCMdict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames) 

        # define additional sub-dictionary to add to PCMdict to house arrays for PCM per class per scheme..
        PCM = {}
        
        # and another to ouse checking tallys
        PCMtally = {}

        # .. for each class within each scheme..
        for cc in xrange(0,Nclasses):
            thisbreakname = breaknames[cc]
            
            # define an empty array for this scheme-class to house PCM
            blankarray = {thisbreakname: cp.deepcopy(zeroMap) }

            # define a zero integer for this scheme-class to house PCMtally
            zerotally = {thisbreakname: 0 }

            # add this blank array to interim PAR dictionary, and initialise the tally counter
            PCM.update(blankarray)
            PCMtally.update(zerotally)

        # add this sub-dictionary to PCMdict for this scheme
        PCM = {'PCM':PCM}
        PCMtally = {'PCMtally':PCMtally}
        PCMdict[scheme].update(PCM)
        PCMdict[scheme].update(PCMtally)

    # loop through all files in 'exportPathDistributed_perpixel' directory....
    dirList = os.listdir(exportPathDistributed_perpixel)
    for i in  xrange(0,len(dirList)):
    
        # does this file contain the '_r' string and is it a gz file - check that this is the extraction output (although should not be anything else in directory)
        fname = dirList[i]
        if fname.find('_r') == -1 | fname.find('.gz') == -1:
            print 'WARNING!! file '+str(fname)+' in '+str(exportPathDistributed_perpixel)+') is not of correct type, trying next file!!!'
            continue

        # if the string looks OK, then import the file
        importarray = np.loadtxt(exportPathDistributed_perpixel+fname)
            
        # check shape of array 
        if importarray.shape != referenceshape:
            print 'WARNING!! file '+str(fname)+ 'has shape '+str(importarray.shape)+' but reference shape is '+str(referenceshape)+', trying next file!!!'
            continue

        # deconstruct filename
        name_parts = deconstructFilename (fname)
        variable = name_parts['variable']

        # if its a meanPR file, add it's running mean values to the global meanPR array
        if variable == str('meanPR'):
            meanPR = meanPR + importarray
            meanPRtally = meanPRtally+1

        # if its a meanPR2 file, add it's running mean values to the global meanPR2 array
        if variable == str('meanPR2'):
            meanPR2 = meanPR2 + importarray
            meanPR2tally = meanPR2tally+1
            
        # if its a meanBUR file, add it's running mean values to the global meanBUR array
        if variable == str('meanBUR'):
            meanBUR = meanBUR + importarray
            meanBURtally = meanBURtally+1

        # if its a meanBUR file, add it's running mean values to the global meanBUR2 array
        if variable == str('meanBUR2'):
            meanBUR2 = meanBUR2 + importarray
            meanBUR2tally = meanBUR2tally+1

        # if its a PCM file
        if variable == str('PCM'):
        
            # get scheme and class ID
            scheme = name_parts['scheme']        
            thisbreakname = name_parts['breakname']

            # add running PCM values to the correct global PCM array for this scheme and class, and update corresponding tally
            PCMdict[scheme]['PCM'][thisbreakname] = PCMdict[scheme]['PCM'][thisbreakname] + importarray
            PCMdict[scheme]['PCMtally'][thisbreakname] = PCMdict[scheme]['PCMtally'][thisbreakname] +1 
            
    # run checks on tallys - they should all be the same for each variable
    if ((meanPRtally == meanPR2tally == meanBURtally == meanBUR2tally)==False):
        print 'WARNING!!! tallys do not match: meanPRtally='+str(meanPRtally)+' meanPR2tally='+str(meanPR2tally)+' meanBURtally='+str(meanBURtally)+' meanBUR2tally='+str(meanBUR2tally)
        
    ## ..loop through each classification scheme 
    for ss in xrange(0,Nschemes): 
        scheme=schemeNames[ss]   
        breaknames = PCMdict[scheme]['BREAKNAMES']
        Nclasses=len(breaknames) 

        # .. for each class within each scheme..
        for cc in xrange(0,Nclasses):
            thisbreakname = breaknames[cc]

            thistally = PCMdict[scheme]['PCMtally'][thisbreakname]
            if(thistally!=meanPRtally):
                print 'WARNING!!! tallys do not match: PCM tally for scheme '+str(scheme)+' class '+str(thisbreakname)+' is '+str(thistally)+' but for meanPRtally is '+str(meanPRtally)
    
    # calculate SD for PR and Burden
    varPR = meanPR2 - np.square(meanPR)
    stdevPR = np.sqrt(varPR)
    varBUR = meanBUR2 - np.square(meanBUR)
    stdevBUR = np.sqrt(varBUR)

    # export mean and SD arrays as asciis
    exportAscii(meanPR,exportPathCombined_perpixel+"meanPR.asc",hdrDict,mask = mask.root.data[:,:])
    exportAscii(stdevPR,exportPathCombined_perpixel+"stdevPR.asc",hdrDict,mask = mask.root.data[:,:])
    exportAscii(meanBUR,exportPathCombined_perpixel+"meanBUR.asc",hdrDict,mask = mask.root.data[:,:])
    exportAscii(stdevBUR,exportPathCombined_perpixel+"stdevBUR.asc",hdrDict,mask = mask.root.data[:,:])

    # for each classification scheme, define an array showing PCM to most likely class (PCMMLC) and what that most likely class is (MLC)
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
            exportAscii(PCMdict[scheme]['PCM'][thisbreakname],exportPathCombined_perpixel+'PCM_'+scheme+'_'+thisbreakname+'.asc',hdrDict,mask = mask.root.data[:,:])

        # export MLC and PCMMLC for this scheme as asciis
        exportAscii(PCMMLC,exportPathCombined_perpixel+'PCMMLC_'+scheme+'.asc',hdrDict,mask = mask.root.data[:,:])
        exportAscii(MLC,exportPathCombined_perpixel+'MLC_'+scheme+'.asc',hdrDict,mask = mask.root.data[:,:])

    return()
#############################################################################################################################################
