# script to download to an instance, before anything executes,  any necessary auxilliary files, and to 
# pre-build any necessary directories

print 'Starting: ECRUNSCRIPT_CONDSIM_PREDOWNLOAD..'

# import libraries
from map_utils import checkAndBuildPaths
from map_utils import S3
from CONDSIM_params import *
import sys

S=S3(keyPath) # initialise key object

# make empty directory on instance to house realisation hdf5 file that wil be generated
print '\n\tBuilding directory: '+realisations_path.rpartition('/')[0]
checkAndBuildPaths(realisations_path.rpartition('/')[0],VERBOSE=True,BUILD=True)

# download from S3 the necessary auxilliary files..

## mcmc trace file
print '\nDownloading burden traces from S3..'
S3bucketname = trace_path.split('/')[-2]
print '\tS3bucketname: '+str(S3bucketname)
S3filename = trace_path.split('/')[-1]
print '\tS3filename: '+str(S3filename)
S.downloadFileFromBucket(S3bucketname,S3filename,trace_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

## global 5km stable mask
print '\nDownloading b5km stable mask from S3..'
S3bucketname = lim5kmbnry_path.split('/')[-2]
print '\tS3bucketname: '+str(S3bucketname)
S3filename = lim5kmbnry_path.split('/')[-1]
print '\tS3filename: '+str(S3filename)
S.downloadFileFromBucket(S3bucketname,S3filename,lim5kmbnry_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

## global 5km urban indicator surface
print '\nDownloading 5km urban indicator surface from S3..'
S3bucketname = urb5km_path.split('/')[-2]
print '\tS3bucketname: '+str(S3bucketname)
S3filename = urb5km_path.split('/')[-1]
print '\tS3filename: '+str(S3filename)
S.downloadFileFromBucket(S3bucketname,S3filename,burdentrace_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

print '\nDONE: ECRUNSCRIPT_extractSummaries_PREDOWNLOAD'
