# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################

import os
import numpy as np
import copy as cp
from rpy import *
from map_utils import quantile_funs as qs
from map_utils import checkAndBuildPaths
from UserDict import UserDict

# import parameters and check filepaths
from extract_params import *
checkAndBuildPaths(filename,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPath,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathCombined,VERBOSE=True,BUILD=True)
checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=True)

# import some r functions
r.source(utilFolder+'GeneralUtility.R')
writeTableWithNamesPY = r['writeTableWithNames']

#############################################################################################################################################
def deconstructFilename (fname):
    # deconstruct filename to assign variables relating to this file 
    part1 = fname.partition('_')
    variable = part1[0]

    if variable==str('PAR'):
       part2 = part1[2].partition('_')
       scheme = part2[0]
       part3 = part2[2].partition('_')
       classname = part3[0]
       part4=part3[2].partition('r')
       part5=part4[2].partition('to')
       startRel = int(part5[0])
       part6=part5[2].partition('.')
       endRel = int(part6[0])
       returnList=variable,scheme,classname,startRel,endRel
       return(returnList)

    if ((variable==str('meanPR')) | (variable==str('BURDEN'))):       
       part2=part1[2].partition('r')
       part3=part2[2].partition('to')
       startRel = int(part3[0])
       part4=part3[2].partition('.')
       endRel = int(part4[0])
       returnList=variable,startRel,endRel
       return(returnList)
       
    print "WARNING!!! file "+fname+" does not contain expected variables"
    return(-9999)

#############################################################################################################################################
def copySubTableToMain(subfname,maintable):

    # read in file
    inputTable = np.loadtxt(exportPath+subfname)

    # define start and end columns on global array    
    name_parts = deconstructFilename(subfname)    
    startRel = name_parts[-2]    
    endRel = name_parts[-1]    
    ncols = inputTable.shape[1]    
    nrel = (endRel-startRel)+1    

    # check dimensions are as expected    
    if inputTable.shape[0] != Nsalb:    
        print( "Warning!!! nrows of "+fname+" is "+str(inputTable.shape[0])+" != expected "+str(Nsalb))    
    if ncols/nrel != n_per:    
        print( "Warning!!! ncols of "+fname+" is "+str(ncols/nrel)+" != expected "+str(n_per))    

    # copy contents of this file to the global array in the relevant position    
    startCol = startRel*n_per    
    endCol = startCol + (nrel*n_per)    
    maintable[:,startCol:endCol:1]=inputTable

    # return modified main table
    return maintable

#############################################################################################################################################
def makeGlobalArray_contVariables(variableName):

    globalarray = cp.deepcopy(blankarray)

    for fname in os.listdir(exportPath):
        if fname.find(variableName+'_')!=-1:

            # copy this array to correct position on global array
            globalarray = copySubTableToMain(fname,globalarray)

    np.savetxt(exportPathCombined+variableName+".txt",globalarray)

    # optionally also create a summary table - mean,SD quantiles etc over all realisations per country
    if summaryStats!=None:
        summObj = getSummariesPerCountry(globalarray)
        summTable = summObj[0]
        summNames = summObj[1]
        writeTableWithNamesPY(summTable,names=summNames,filepath=exportPathCombined+variableName+"_summary.csv",SEP=",")

    return(globalarray) 

#############################################################################################################################################
def makeGlobalArray_categoVariables(variableName):

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
            for fname in os.listdir(exportPath):
                if (fname.find(variableName+'_')!=-1) & (fname.find(scheme)!=-1) & (fname.find(thisbreakname)!=-1) :

                    # copy this array to correct position on global array
                    globalarray = copySubTableToMain(fname,globalarray)

            np.savetxt(exportPathCombined+variableName+'_'+scheme+'_'+thisbreakname+'.txt',globalarray)

            # optionally also create a summary table - mean,SD quantiles etc over all realisations per country
            if summaryStats!=None:
                summObj = getSummariesPerCountry(globalarray)
                summTable = summObj[0]
                summNames = summObj[1]
                writeTableWithNamesPY(summTable,names=summNames,filepath=exportPathCombined+variableName+'_'+scheme+'_'+thisbreakname+"_summary.csv",SEP=",")

#############################################################################################################################################
def getSummariesPerCountry(globalarray):

    # read in array of salb IDs
    uniqueSalb=fromfile(uniqueSalb_path,sep=",")
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


## loop through all files in 'exportPath' directory, get list of unique realization ranges, and thereby
## calculate total top-level realizations present accross output files in directory

relRanges=list()  
i=0
for fname in os.listdir(exportPath):

    # does this file contain the '_r' string and is it a text file - check that this is the extraction output (although should not be anything else in directory)
    if fname.find('_r') == -1 | fname.find('.txt') == -1: continue

    # record name of first file for later use
    if i==0:
        firstfname = fname
        i=1

    # deconstruct filename
    name_parts = deconstructFilename (fname)
    relRange=name_parts[-2],name_parts[-1]

    relRanges.append(relRange)

# get array of unique rel ranges
uniqueRelRanges=np.unique(relRanges)
Nunique = uniqueRelRanges.shape[0] 

# calculate total number of top-level realisations present accross these files
n_realizations_infiles=0
for i in xrange(0,Nunique):
    subtotal = (uniqueRelRanges[i][1] - uniqueRelRanges[i][0])+1
    n_realizations_infiles = n_realizations_infiles + subtotal

# import first file and obtain dimensions to reveal N_uniqueSalb and n_per 
name_parts = deconstructFilename (firstfname)
tableshape = np.loadtxt(exportPath+firstfname).shape
Nsalb = tableshape[0]
ncols = tableshape[1] 
n_per = ncols/((name_parts[-1]-name_parts[-2])+1)

# define generic blank array to house combined tables for each variable
blankarray = np.repeat(-9999.,n_realizations_infiles*n_per*Nsalb).reshape(Nsalb,n_per*n_realizations_infiles)

# loop through all meanPR_* files and add them to global array, then export global array, and optionallly accopmanying summary table
temp=makeGlobalArray_contVariables('meanPR')

# loop through all BURDEN_* files and add them to global array, then export global array, and optionallly accopmanying summary table
temp=makeGlobalArray_contVariables('BURDEN')

# loop through all PAR_* files and add them to global arrays, then export global arrays, and optionallly accopmanying summary table
makeGlobalArray_categoVariables('PAR')
