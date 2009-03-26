import numpy as np
import tables as tb
import copy as cp
from rpy import *
from extract_PYlib import examineSalb

# import some r functions
from extract_params import utilFolder
r.source(utilFolder+'GeneralUtility.R')
writeTableWithNamesPY = r['writeTableWithNames']

# set some parameters
salb1km_path = "/home/pwg/mbg-world/datafiles/auxiliary_data/salb1km-e_AF.hdf5"
gr001km_path= "/home/pwg/mbg-world/datafiles/auxiliary_data/gr001km_AF.hdf5"
lims1km_path= "/home/pwg/mbg-world/datafiles/auxiliary_data/lims1km-e_AF.hdf5"
outputTable_path = "/home/pwg/mbg-world/extraction/FixedPopulationExtraction.csv"

# check paths
from map_utils import checkAndBuildPaths
checkAndBuildPaths(salb1km_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(lims1km_path,VERBOSE=True,BUILD=False)


# open link to salb grid, 3-level limits grid, and population grid
salb1km = tb.openFile(salb1km_path, mode = "r")    
gr001km = tb.openFile(gr001km_path, mode = "r")    
lims1km = tb.openFile(lims1km_path, mode = "r") 

# run check that input grids are the same shape
if(((np.shape(salb1km.root.data)==np.shape(gr001km.root.data)==np.shape(lims1km.root.data))==False)):
    print "WARNING!! input grids are of uneven shape. salb1km="+str(np.shape(salb1km.root.data))+"; gr001km="+str(np.shape(gr001km.root.data))+"; lims1km="+str(np.shape(lims1km.root.data))

# run extract salb to get list of unique salb IDs and corresponding pixel count
salbDict = examineSalb (salb1km)
uniqueSalb = salbDict['uniqueSalb']
countSalb = salbDict['count']
Nsalb=len(uniqueSalb) 

# intialise empty arrays to house total population in each country in each limits class
Pop_all = zeros(Nsalb) 
Pop_nullrisk = zeros(Nsalb)
Pop_riskfree = zeros(Nsalb)
Pop_unstable = zeros(Nsalb)
Pop_stable = zeros(Nsalb)
countSalb_check = zeros(Nsalb) #  to house counted pixels to cross-check with 'countSalb' returned from examineSalb

interimCnt = 0

n_rows = len(salb1km.root.lat)
for jj in xrange(0,n_rows):

    interimCnt = interimCnt+1
    if interimCnt==100:
        print "On row "+str(jj)+" of "+str(n_rows)
        interimCnt=0

    # get this row of 1km Salb,population surface,and limits (assumes they are correct way up i.e. map view)  
    salb1km_ROW = salb1km.root.data[jj,:]
    gr001km_ROW = gr001km.root.data[jj,:]
    lims1km_ROW = lims1km.root.data[jj,:]

    # define a blank array of zeroes of same size as 1km row - that will be duplicated for various uses later
    zeroChunk = zeros(product(gr001km_ROW.shape)).reshape(gr001km_ROW.shape)

    # define ID matrices for this row for each limits class
    nullriskIDmatrix = cp.deepcopy(zeroChunk)
    riskfreeIDmatrix = cp.deepcopy(zeroChunk)
    unstableIDmatrix = cp.deepcopy(zeroChunk)
    stableIDmatrix = cp.deepcopy(zeroChunk)

    nullriskIDmatrix[lims1km_ROW==-9999]=1    
    riskfreeIDmatrix[lims1km_ROW==0]=1    
    unstableIDmatrix[lims1km_ROW==1]=1
    stableIDmatrix[lims1km_ROW==2]=1

    # how many unique salb IDs in these rows (after removing -9999 cells from sea and/or non-stable areas)
    uniqueSalb_ROW = unique(salb1km_ROW)
    Nsalb_ROW =len(uniqueSalb_ROW)

    # which rows on the all-country arrays correpsond to these salb IDs?
    salbLUT_ROW=array(r.match(uniqueSalb_ROW,uniqueSalb))
    #.. take care of special case of only 1 unique country ID in this slice - must ensure salbLUT_ROW is an array even if onluy one element
    if Nsalb_ROW==1:
         salbLUT_ROW.shape=1
    salbLUT_ROW = salbLUT_ROW-1  # converts from r to python naming convention   

    # loop through each country present in this chunk
    for rr in xrange(0,Nsalb_ROW): 

        # which row on the all-country arrays does this country correpsond to?        
        whichRow = salbLUT_ROW[rr]

        # which pixels are in this country
        countryIDmatrix = cp.deepcopy(zeroChunk)
        countryIDmatrix[salb1km_ROW==uniqueSalb_ROW[rr]]=1

        # add to count of country pixels (for later cross-checking)
        countSalb_check[whichRow] = countSalb_check[whichRow] + sum(countryIDmatrix)

        # calculate sum of population in each limits class for this country/row
        Pop_all_ROW = sum(countryIDmatrix * gr001km_ROW)
        Pop_nullrisk_ROW = sum(countryIDmatrix * nullriskIDmatrix * gr001km_ROW)
        Pop_riskfree_ROW = sum(countryIDmatrix * riskfreeIDmatrix * gr001km_ROW)
        Pop_unstable_ROW = sum(countryIDmatrix * unstableIDmatrix * gr001km_ROW)
        Pop_stable_ROW = sum(countryIDmatrix * stableIDmatrix * gr001km_ROW)

        # add these sums to global population count vectors for this country:
        Pop_all[whichRow] = Pop_all[whichRow] + Pop_all_ROW
        Pop_nullrisk[whichRow] = Pop_nullrisk[whichRow] + Pop_nullrisk_ROW
        Pop_riskfree[whichRow] = Pop_riskfree[whichRow] + Pop_riskfree_ROW
        Pop_unstable[whichRow] = Pop_unstable[whichRow] + Pop_unstable_ROW
        Pop_stable[whichRow] = Pop_stable[whichRow] + Pop_stable_ROW

        ## run checks on these sums:
        # should be no 'null' risk pixels in areas of non-null salb ID
        if (uniqueSalb_ROW[rr] != -9999):
            if Pop_nullrisk_ROW > 0:
                print "WARNING!! "+str(sum(countryIDmatrix * nullriskIDmatrix))+" pixels and "+str(Pop_nullrisk_ROW)+" population found in null risk area in salbID "+str(uniqueSalb_ROW[rr])+" , in row "+str(jj)  

        # should be no non-zero population pixels in areas of null salb ID
        if (uniqueSalb_ROW[rr] == -9999):
            if sum(countryIDmatrix * gr001km_ROW) > 0:
                print "WARNING!! "+str(sum(countryIDmatrix * gr001km_ROW))+" population found in null salb (-9999) area in row "+str(jj)  

# run various checks on final output for each country
for ii in xrange(0,Nsalb):

    # sum of pixels in each country matches input countSalb vector
    if countSalb_check[ii]!=countSalb[ii]: print "WARNING!! sum of pixels is inconsistent (countSalb_check="+str(countSalb_check[ii])+" countSalb="+str(countSalb[ii])+") in country salbID = "+str(uniqueSalb[ii])

    # sum of population in each limits class matches total population in each country 
    popsumclasses = Pop_nullrisk[ii] + Pop_riskfree[ii] + Pop_unstable[ii] + Pop_stable[ii]
    popsum = Pop_all[ii]
    if (popsumclasses!=popsum):
        print "WARNING!! sum of population accross classes ("+str(popsumclasses)+") != total population ("+str(popsum)+") in country salbID = "+str(uniqueSalb[ii])

# build and export array of country population summaries
pop_table = zeros(Nsalb*6).reshape(Nsalb,6)
pop_table[:,0] = uniqueSalb
pop_table[:,1] = Pop_all
pop_table[:,2] = Pop_nullrisk
pop_table[:,3] = Pop_riskfree
pop_table[:,4] = Pop_unstable
pop_table[:,5] = Pop_stable
pop_table_names = list(['salbID','Pop_all','Pop_nullrisk','Pop_riskfree','Pop_unstable','Pop_stable'])
writeTableWithNamesPY(pop_table,names=pop_table_names,filepath=outputTable_path,SEP=",")
