# make sure we are in the correct directory
cd /root/mbg-world/mbgw-scripts/

# import libraries
from map_utils import checkAndBuildPaths
from extract_PYlib import *
from boto_PYlib import *


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

BURDEN = str(sys.argv[7])
if BURDEN == 'None': BURDEN = None
else: BURDEN = int(BURDEN)

# build realisation block import path
hdf5block = realisations_path
hdf5block = hdf5block.replace('FILESTARTREL',str(FileStartRel))
hdf5block = hdf5block.replace('FILEENDREL',str(FileEndRel))

# download this realisation file from S3 storage
S3bucketname = hdf5block.split('/')[-2]
S3filename = hdf5block.split('/')[-1]
downloadFileFromBucket(S3bucketname,S3filename,hdf5block,overwriteContent=False,makeDirectory=True,VERBOSE=True)
checkAndBuildPaths(hdf5block,VERBOSE=False,BUILD=False)

# download from S3 the other necessary file (optionally need 5km grump for burden map)
if (BURDEN==True):
    S3bucketname = gr005km_path.split('/')[-2]
    S3filename = gr005km_path.split('/')[-1]
    downloadFileFromBucket(S3bucketname,S3filename,gr005km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
    checkAndBuildPaths(gr005km_path,VERBOSE=False,BUILD=False)

# check path for exports exists
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)

# now call extractSummaries_perpixel substituting in the formatted sys args 
extractSummaries_perpixel (slices,a_lo,a_hi,n_per,FileStartRel,FileEndRel,totalN,startRel=None,endRel=None,BURDEN=False)