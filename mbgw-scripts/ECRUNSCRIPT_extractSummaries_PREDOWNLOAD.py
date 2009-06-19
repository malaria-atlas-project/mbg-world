# script to download to an instance, before anything execute,  any necessary auxilliary files, and to 
# pre-build any necessary directories

print 'Starting: ECRUNSCRIPT_extractSummaries_PREDOWNLOAD..'

# import libraries
from map_utils import checkAndBuildPaths
from  map_utils import boto_PYlib
from extract_params import *
import sys

S=S3(keyPath) # initialise key object

# deal with system arguments
BURDEN = True
PERPIXEL = True
PERCOUNTRY = True

if sys.argv[1] == 'False' : BURDEN=False
if sys.argv[2] == 'False' : PERPIXEL=False
if sys.argv[3] == 'False' : PERCOUNTRY=False

# make empty directory on instance to house realisation hdf5 file downloaded from S3
print '\n\tBuilding directory: '+realisations_path.rpartition('/')[0]
checkAndBuildPaths(realisations_path.rpartition('/')[0],VERBOSE=True,BUILD=True)


# optionally download the burden traces from S3 storage
if (BURDEN==True):
    print '\nDownloading burden traces from S3..'
    S3bucketname = burdentrace_path.split('/')[-2]
    print '\tS3bucketname: '+str(S3bucketname)
    S3filename = burdentrace_path.split('/')[-1]
    print '\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,burdentrace_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

if (PERPIXEL==True):
    # make empty directory on instance to house output files ready to be uploaded back to S3
    print '\n\tBuilding directory: '+exportPathDistributed_perpixel
    checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)    

    if (BURDEN==True):
        print '\n\tDownloading auxilliary file to '+str(grump5km_path)+' from S3..'
        S3bucketname = grump5km_path.split('/')[-2]
        print '\t\tS3bucketname: '+str(S3bucketname)
        S3filename = grump5km_path.split('/')[-1]
        print '\t\tS3filename: '+str(S3filename)
        S.downloadFileFromBucket(S3bucketname,S3filename,grump5km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

if (PERCOUNTRY==True):
    # make empty directory on instance to house output files ready to be uploaded back to S3
    print '\n\tBuilding directory: '+exportPathDistributed_country
    checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=True)  

    print '\n\tDownloading auxilliary file to '+str(grump1km_path)+' from S3..'
    S3bucketname = grump1km_path.split('/')[-2]
    print '\t\tS3bucketname: '+str(S3bucketname)
    S3filename = grump1km_path.split('/')[-1]
    print '\t\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,grump1km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

    print '\n\tDownloading auxilliary file to '+str(salblim1km_path)+' from S3..'    
    S3bucketname = salblim1km_path.split('/')[-2]    
    print '\t\tS3bucketname: '+str(S3bucketname)    
    S3filename = salblim1km_path.split('/')[-1]
    print '\t\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,salblim1km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)

print '\nDONE: ECRUNSCRIPT_extractSummaries_PREDOWNLOAD'
