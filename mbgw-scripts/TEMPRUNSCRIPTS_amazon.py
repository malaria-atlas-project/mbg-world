
#################################################################################
EXTRACT PER-COUNTRY MEAN PR,BURDEN,PAR

from extract_PYlib import *

# check filepaths stated in parameter file
from map_utils import checkAndBuildPaths
#checkAndBuildPaths(filename,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathCombined_country,VERBOSE=True,BUILD=True)
checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=True)

#a=time.time()
#extractSummaries_country([slice(None,None,None), slice(None,None,None), slice(0,12,None)],2,10,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
#print "all done from PYlib"
#print("TOTAL TIME: "+(str(time.time()-a)))
#OR

extractSummaries_country([slice(None,None,None), slice(None,None,None), slice(0,12,None)],2,10,1,1,2)


#################################################################################
EXTRACT PER-PIXEL PR, CLASS, AND BURDEN SUMMARIES
#################################################################################
COMBINE DISTRIBUTED COUNTRY AND PER PIXEL EXTRACTIONS

from extract_combineExtractions import *

from map_utils import checkAndBuildPaths
#checkAndBuildPaths(filename,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=False)
checkAndBuildPaths(exportPathCombined_country,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=False)
checkAndBuildPaths(exportPathCombined_perpixel,VERBOSE=True,BUILD=True)
checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(lim5kmbnry_path,VERBOSE=True,BUILD=False)

#temp=combineDistribExtractions_perpixel()
combineDistribExtractions_country()
#################################################################################


from extract_params import *
from map_utils import checkAndBuildPaths


checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathCombined_country,VERBOSE=True,BUILD=True)
checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=False)
checkAndBuildPaths(exportPathCombined_country,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=False)
checkAndBuildPaths(exportPathCombined_perpixel,VERBOSE=True,BUILD=True)
checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(gr001km_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(lim5kmbnry_path,VERBOSE=True,BUILD=False)
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)
checkAndBuildPaths(exportPathCombined_perpixel,VERBOSE=True,BUILD=True)
checkAndBuildPaths(lim5kmbnry_path,VERBOSE=True,BUILD=True)
checkAndBuildPaths(STDOUTPUT,VERBOSE=True,BUILD=True)














