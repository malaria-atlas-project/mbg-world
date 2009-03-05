# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################

import os
import numpy as np
from extract_params import *
import copy as cp

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

    if variable==str('meanPR'):       
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
    name_parts = deconstructFilename (subfname)    
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

    np.savetxt(exportPath+variableName+".txt",globalarray) 

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

            np.savetxt(exportPath+variableName+'_'+scheme+'_'+thisbreakname+'.txt',globalarray)

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

# loop through all meanPR_* files and add them to global array, then export global array
makeGlobalArray_contVariables('meanPR')

# loop through all PAR_* files and add them to global arrays, then export global arrays
makeGlobalArray_categoVariables('PAR')
