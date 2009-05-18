# script to download to an instance, before anything execute,  any necessary auxilliary files, and to 
# pre-build any necessary directories

print 'Starting: ECRUNSCRIPT_extractSummaries_perpixel_PREDOWNLOAD..'

# import libraries
from map_utils import checkAndBuildPaths
from boto_PYlib import *
import sys

S=S3() # initialise key object

# deal with system arguments
REGION = str(sys.argv[1])
BURDEN = bool(sys.argv[2])
PERPIXEL = bool(sys.argv[3])
PERCOUNTRY = bool(sys.argv[4])

if REGION=="AF": from extract_params_AF import * 
if REGION=="AM": from extract_params_AM import *
if REGION=="AS": from extract_params_AS import *

# make empty directory on instance to house realisation hdf5 file downloaded from S3
print '\n\tBuilding directory: '+realisations_path.rpartition('/')[0]
checkAndBuildPaths(realisations_path.rpartition('/')[0],VERBOSE=True,BUILD=True)

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
