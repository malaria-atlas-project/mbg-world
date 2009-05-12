# check in correct dirctory

# import libraries
from map_utils import checkAndBuildPaths
from extract_PYlib import *
from boto_PYlib import *
from extract_params import *
import os

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
else: BURDEN = bool(BURDEN)

# build realisation block import path
hdf5block_path = realisations_path
hdf5block_path = hdf5block_path.replace('FILESTARTREL',str(FileStartRel))
hdf5block_path = hdf5block_path.replace('FILEENDREL',str(FileEndRel))

print "from RUNSCRIPT_extractSummaries_perpixel:\n"
print "n_per: "+str(n_per)
print "FileStartRel: "+str(FileStartRel)
print "FileEndRel: "+str(FileEndRel)
print "totalN: "+str(totalN)
print "startRel: "+str(startRel)
print "endRel: "+str(endRel)
print "BURDEN: "+str(BURDEN)

# download this realisation file from S3 storage
print '\nDownloading realisation from S3..'
S3bucketname = hdf5block_path.split('/')[-2]
print '\tS3bucketname: '+str(S3bucketname)
S3filename = hdf5block_path.split('/')[-1]
print '\tS3filename: '+str(S3filename)
S.downloadFileFromBucket(S3bucketname,S3filename,hdf5block_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
checkAndBuildPaths(hdf5block_path,VERBOSE=False,BUILD=False)

# download from S3 the other necessary file (optionally need 5km grump for burden map)
if (BURDEN==True):
    print '\nDownloading other files from S3..'
    S3bucketname = gr005km_path.split('/')[-2]
    print '\tS3bucketname: '+str(S3bucketname)
    S3filename = gr005km_path.split('/')[-1]
    print '\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,gr005km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
    checkAndBuildPaths(gr005km_path,VERBOSE=False,BUILD=False)

# check path for exports exists
print '\nchecking path for export exists..'
checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)

# now call extractSummaries_perpixel substituting in the formatted sys args 
print '\nrunning extractSummaries_perpixel..'
extractSummaries_perpixel ([slice(None,None,None), slice(None,None,None), slice(0,12,None)],2,10,n_per,FileStartRel,FileEndRel,totalN,startRel,endRel,BURDEN)

# now upload the output back to the S3 storage
print '\nrunning extractSummaries_perpixel..'
#S.uploadDirectoryAsBucket('distributedoutput_perpixel',exportPathDistributed_perpixel,uploadConstituentFiles=True,overwriteContent=True)

## loop through all files in local export storage
for fname in os.listdir(exportPathDistributed_country):

    # if file contains string indicating an outcome from this realisation set, then upload to S3
    if (fname.find('r'+str(FileStartRel)+'to'+str(FileEndRel)+'.gz')>0):
        filepath = exportPathDistributed_country+fname
        print 'uploading '+filepath+' to S3 bucket '+'distributedoutput_perpixel'
        S.uploadFileToBucket('distributedoutput_perpixel',filePath,overwriteContent=True,makeBucket=True,VERBOSE=True)
