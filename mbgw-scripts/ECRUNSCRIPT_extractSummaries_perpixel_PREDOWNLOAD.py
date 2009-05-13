print 'Starting: ECRUNSCRIPT_extractSummaries_perpixel_PREDOWNLOAD..'

# import libraries
from map_utils import checkAndBuildPaths
from boto_PYlib import *
from extract_params import *

S=S3() # initialise key object

# deal with system arguments
BURDEN = str(sys.argv[1])
if BURDEN == 'None': BURDEN = None  
else: BURDEN = bool(BURDEN)

# download from S3 the other necessary file (optionally need 5km grump for burden map)
if (BURDEN==True):
    print '\n\tDownloading other files from S3..'
    S3bucketname = grump5km_path.split('/')[-2]
    print '\t\tS3bucketname: '+str(S3bucketname)
    S3filename = grump5km_path.split('/')[-1]
    print '\t\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,grump5km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
    checkAndBuildPaths(grump5km_path,VERBOSE=False,BUILD=False)

# make empty directory on instance to house realisation hdf5 file downloaded from S3
print '\n\tBuilding directory: '+realisations_path.split('/')[-2]
checkAndBuildPaths(realisations_path.split('/')[-2],VERBOSE=True,BUILD=True)

# make empty directory on instance to house output files ready to be uploaded back to S3
print '\n\tBuilding directory: '+exportPathDistributed_perpixel
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)

print '\nDONE: ECRUNSCRIPT_extractSummaries_perpixel_PREDOWNLOAD'
