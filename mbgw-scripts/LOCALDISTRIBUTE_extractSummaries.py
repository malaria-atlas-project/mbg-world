# example command line:
# python LOCALDISTRIBUTE_extractSummaries extract_params_AF.py 

# import libraries
import numpy as np
from map_utils import checkAndBuildPaths
import time
import sys
import shutil
import os


# deal with system argument defining version of parameter file - create copy with generic name extract_params.py
PARAMFILE = sys.argv[1]  ## defines name of python file housing the parmeter definitions (e.g. extract_params_AF.py)
print 'defining local paramfile "extract_params.py" from '+str(PARAMFILE.partition('.')[0])
shutil.copy(PARAMFILE,'extract_params.py')
from extract_params import *
from extract_PYlib import *

# set path to realisation folder an extract folder path and  generic file name
folderPath = realisations_path.rpartition('/')[0]
relPath = realisations_path.rsplit('/')[-1]

# call queryRealizationsInFolder to obtain number and start/end realisation numbers of these realisation files
relDict = queryRealizationsInFolder(folderPath,relPath,VERBOSE=True)
print '\nquerying folder '+str(folderPath)+' : found '+str(relDict['Nrealisations'])+' realisations accross '+str(relDict['Nfiles'])+' files.'
if relDict['Nrealisations']==0: raise RuntimeError ('!!No realisations to extract so quitting!!')

# set realization number parameters
NRELS = relDict['Nrealisations']
NJOBS = relDict['Nfiles']

####################################TEMP
#NJOBS = 1
#NRELS = 1
####################################TEMP

FileStartRels = relDict['StartRelList']
FileEndRels = relDict['EndRelList']
NPER  = 1
NTOTALREL = NRELS*NPER

####################################TEMP 
#NTOTALREL = 1
####################################TEMP

# call examineSalb to provide lists of unique spatial units, and optionally their areas and populations
temp=examineSalb(salb=salb5km_path, grump=grump5km_path,salbpop_path=salbpop_path, pixarea=pixarea5km_path,salbarea_path=salbarea_path,uniqueSalb_path=uniqueSalb_path,pixelN_path=pixelN_path,ignore=np.array([-9999]))

## main job

for i in xrange(NJOBS):

    print 'Running extractions for realisation '+str(i)+' of '+str(NJOBS)

    # build filename of hdf5 realization file
    hdf5block_path = realisations_path
    hdf5block_path = hdf5block_path.replace('FILESTARTREL',str(FileStartRels[i]))
    hdf5block_path = hdf5block_path.replace('FILEENDREL',str(FileEndRels[i]))
    
    if PERPIXEL is True:

        print '\n running PERPIXEL extraction:'

        # check path for per-pixel exports exists
        print '\n\nchecking path for export exists..'
        checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)

        # now call extractSummaries_perpixel substituting in the formatted sys args 
        print '\n\nrunning extractSummaries_perpixel..'
        extractSummaries_perpixel ([slice(None,None,None), slice(None,None,None), MonthsSlice],a_lo,a_hi,NPER,FileStartRels[i],FileEndRels[i],NTOTALREL,None,None,do_PRMap,do_BurdenMap,do_RoMap)


    if PERCOUNTRY is True:

        print '\n running PERCOUNTRY extraction:'

        # check path for per-country exports exists
        print '\nchecking path for export exists..'
        checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=True)

        # now call extractSummaries_country substituting in the formatted sys args 
        print '\nrunning extractSummaries_country..'
        extractSummaries_country([slice(None,None,None), slice(None,None,None), MonthsSlice],a_lo,a_hi,NPER,FileStartRels[i],FileEndRels[i],None,None, do_AREALMEANPR,do_POPMEANPR,do_PAR,do_BURDEN,do_AREALMEANRo,do_POPMEANRo)



