## example call:
# run LOCALDISTRIBUTE_combineDistribExtractions extract_params_AM.py

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
from extract_combineExtractions import *

if PERPIXEL is True:

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

    # build path for output to house combined per-pixel output maps
    checkAndBuildPaths(exportPathCombined_perpixel,VERBOSE=True,BUILD=True)

    checkAndBuildPaths(limbnry_path,VERBOSE=True,BUILD=False)

    if (do_BurdenMap==True): checkAndBuildPaths(grump_path,VERBOSE=True,BUILD=False)

    # now call extractSummaries_perpixel substituting in the formatted sys args 
    print '\n\tCalling combineDistribExtractions_perpixel'
    combineDistribExtractions_perpixel()

    # now upload the output back to the S3 storage

if PERCOUNTRY is True:

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

    checkAndBuildPaths(salblim_path,VERBOSE=True,BUILD=False)
    checkAndBuildPaths(salb_path,VERBOSE=True,BUILD=False)
    
    # build paths to directory to house uniqueSalb.txt, pixelN.txt,uniqueSalbAllCountries.txt, pixelNAllCountries.txt 
    checkAndBuildPaths(uniqueSalb_path.rpartition('/')[0],VERBOSE=True,BUILD=True)
    
    # check paths for examinSalb output
    checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=False)
    checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=False)
    if (do_BURDEN | do_POPMEANPR | do_POPMEANRo | do_PAR): checkAndBuildPaths(salbpop_path,VERBOSE=True,BUILD=False)
    if (do_AREALMEANPR | do_AREALMEANRo):checkAndBuildPaths(salbarea_path,VERBOSE=True,BUILD=False)

    # build path for output to house combined per-pixel output maps
    checkAndBuildPaths(exportPathCombined_country,VERBOSE=True,BUILD=True)

    # now call extractSummaries_country substituting in the formatted sys args 
    print '\n\tCalling combineDistribExtractions_country'
    combineDistribExtractions_country()
            

