# check in correct dirctory

# import libraries
from map_utils import checkAndBuildPaths
from extract_combineDistribExtractions import *
from boto_PYlib import *
from extract_params import *

S=S3() # initialise key object

# deal with system arguments

############################TEMP
#n_per = 3
#FileStartRel = 0
#FileEndRel = 1
#totalN = 3
#startRel = None
#endRel = None
#BURDEN = True
################################

call: combineDistribExtractions_perpixel()

print "from ECRUNSCRIPT_combineDistribExtractions_perPixel:\n"

# download from S3 contents of bucket 'distributedoutput_perpixel', will automatically build the local directory if necessary
S.downloadBucketContents('distributedoutput_perpixel',exportPathDistributed_perpixel,overwriteContent=FALSE,VERBOSE=True)

# build path for output to house combined per-pixel output maps
checkAndBuildPaths(exportPathCombined_perpixel,VERBOSE=TRUE,BUILD=TRUE)

# download from S3 the other necessary files (optionally need 5km grump for burden map)
S3bucketname = lim5kmbnry_path.split('/')[-2]
S3filename = lim5kmbnry_path.split('/')[-1]
S.downloadFileFromBucket(S3bucketname,S3filename,lim5kmbnry_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
checkAndBuildPaths(lim5kmbnry_path,VERBOSE=False,BUILD=False)

# check path for exports exists
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)

# now call extractSummaries_perpixel substituting in the formatted sys args 
combineDistribExtractions_perpixel()

# now upload the output back to the S3 storage
S.uploadDirectoryAsBucket('CombinedOutput_perpixel',exportPathCombined_perpixel,uploadConstituentFiles=True,overwriteContent=True)
