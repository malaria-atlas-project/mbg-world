
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

# make sure we are in the correct directory
cd /root/mbg-world/mbgw-scripts/

# deal with system arguments
n_per = int(sys.argv[1])
FileStartRel = int(sys.argv[2])
FileEndRel = int(sys.argv[3])
totalN = int(sys.argv[4])

startRel = str(sys.argv[5])
if startRel == 'None': startRel = None
else: startRel = int(startRel)

endRel = str(sys.argv[6])
if endRel == 'None': endRel = None
else: endRel = int(endRel)

# import extraction function library 
from extract_PYlib import *

# build realisation block import path
filename = realisations_path
filename = filename.replace('FILESTARTREL',str(FileStartRel))
filename = filename.replace('FILEENDREL',str(FileEndRel))

# download this realisation file from S3 storage
from boto_PYlib import *
S3bucketname = filename.split('/')[-2]
S3filename = filename.split('/')[-1]
downloadFileFromBucket(S3bucketname,S3filename,filename,overwriteContent=True,makeDirectory=True,VERBOSE=True)



checkAndBuildPaths(filename,VERBOSE=False,BUILD=False)

# check relevant filepaths stated in parameter file (although the hdf5 block path is tested within the function)
from map_utils import checkAndBuildPaths
checkAndBuildPaths(gr005km_path,VERBOSE=True,BUILD=FALSE)
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)


slices,a_lo,a_hi,n_per,FileStartRel,FileEndRel,totalN,startRel=None,endRel=None,BURDEN=False


extractSummaries_perpixel ([slice(None,None,None), slice(None,None,None), slice(0,12,None)],2,10,1,0,1,totalN=1,BURDEN=True)

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














