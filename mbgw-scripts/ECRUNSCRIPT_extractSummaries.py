print "from RUNSCRIPT_extractSummaries:\n"

# import libraries
from map_utils import checkAndBuildPaths
from extract_PYlib import *
#from boto_PYlib import *
from map_utils import S3
import os
from extract_params import *
import sys

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

BURDEN = True
PERPIXEL = True
PERCOUNTRY = True

if sys.argv[7] == 'False' : BURDEN=False
if sys.argv[8] == 'False' : PERPIXEL=False
if sys.argv[9] == 'False' : PERCOUNTRY=False

S=S3(keyPath) # initialise key object


# build realisation block import path
hdf5block_path = realisations_path
hdf5block_path = hdf5block_path.replace('FILESTARTREL',str(FileStartRel))
hdf5block_path = hdf5block_path.replace('FILEENDREL',str(FileEndRel))

print "n_per: "+str(n_per)
print "FileStartRel: "+str(FileStartRel)
print "FileEndRel: "+str(FileEndRel)
print "totalN: "+str(totalN)
print "startRel: "+str(startRel)
print "endRel: "+str(endRel)
print "BURDEN: "+str(BURDEN) 
print "PERPIXEL: "+str(PERPIXEL) 
print "PERCOUNTRY: "+str(PERCOUNTRY) 

# download this realisation file from S3 storage
print '\nDownloading realisation from S3..'
S3bucketname = hdf5block_path.split('/')[-2]
print '\tS3bucketname: '+str(S3bucketname)
S3filename = hdf5block_path.split('/')[-1]
print '\tS3filename: '+str(S3filename)
S.downloadFileFromBucket(S3bucketname,S3filename,hdf5block_path,overwriteContent=False,makeDirectory=False,VERBOSE=True)

if PERPIXEL==True:

    print '\n running PERPIXEL extraction'

    # check path for per-pixel exports exists
    print '\nchecking path for export exists..'
    checkAndBuildPaths(exportPathDistributed_perpixel,VERBOSE=True,BUILD=True)

    # now call extractSummaries_perpixel substituting in the formatted sys args 
    print '\nrunning extractSummaries_perpixel..'
    extractSummaries_perpixel ([slice(None,None,None), slice(None,None,None), MonthsSlice],2,10,n_per,FileStartRel,FileEndRel,totalN,startRel,endRel,BURDEN)

    # now upload the output back to the S3 storage
    #S.uploadDirectoryAsBucket('distributedoutput_perpixel',exportPathDistributed_perpixel,uploadConstituentFiles=True,overwriteContent=True)

    ## loop through all files in local export storage
    for fname in os.listdir(exportPathDistributed_perpixel):

        # if file contains string indicating an outcome from this realisation set, then upload to S3 (try a few times before failing this section)
        #if (fname.find('r'+str(FileStartRel)+'to'+str(FileEndRel)+'.gz')>0):
        if True:
            filepath = exportPathDistributed_perpixel+fname
            print 'uploading '+filepath+' to S3 bucket '+ exportPathDistributed_perpixel.split('/')[-2].lower()
            failCount = 0
            while failCount<=3:
                try:
                    S.uploadFileToBucket(exportPathDistributed_perpixel.split('/')[-2].lower(),filepath,overwriteContent=True,makeBucket=True,VERBOSE=True)
                    break
                except RuntimeError:
                    failCount+=1 
                    if failCount<=3:
                        print 'upload of '+filepath+' to S3 bucket '+exportPathDistributed_perpixel.split('/')[-2].lower()+' failed '+str(failCount) +' times: retrying..'
                    else:
                        print 'upload of '+filepath+' to S3 bucket '+exportPathDistributed_perpixel.split('/')[-2].lower()+' failed '+str(failCount) +' times: GIVING UP - FILE WILL NOT BE COPIED!!'  


if PERCOUNTRY==True:

    print '\n running PERCOUNTRY extraction'
                    
    # check path for per-country exports exists
    print '\nchecking path for export exists..'
    checkAndBuildPaths(exportPathDistributed_country,VERBOSE=True,BUILD=True)

    # now call extractSummaries_country substituting in the formatted sys args 
    print '\nrunning extractSummaries_country..'
    #extractSummaries_country([slice(None,None,None), slice(None,None,None), slice(252,264,None)],2,10,n_per,FileStartRel,FileEndRel,startRel,endRel)
    extractSummaries_country([slice(None,None,None), slice(None,None,None), MonthsSlice],2,10,n_per,FileStartRel,FileEndRel,startRel,endRel)

    ## loop through all files in local export storage
    for fname in os.listdir(exportPathDistributed_country):

        # if file contains string indicating an outcome from this realisation set, then upload to S3 (try a few times before failing this section)
        #if (fname.find('r'+str(FileStartRel)+'to'+str(FileEndRel)+'.txt')>0):
        if True:
            filepath = exportPathDistributed_country+fname
            print 'uploading '+filepath+' to S3 bucket '+ exportPathDistributed_country.split('/')[-2].lower()
            failCount = 0
            while failCount<=3:
                try:
                    S.uploadFileToBucket(exportPathDistributed_country.split('/')[-2].lower(),filepath,overwriteContent=True,makeBucket=True,VERBOSE=True)
                    break
                except RuntimeError:
                    failCount+=1 
                    if failCount<=3:
                        print 'upload of '+filepath+' to S3 bucket '+exportPathDistributed_country.split('/')[-2].lower()+' failed '+str(failCount) +' times: retrying..'
                    else:
                        print 'upload of '+filepath+' to S3 bucket '+exportPathDistributed_country.split('/')[-2].lower()+' failed '+str(failCount) +' times: GIVING UP - FILE WILL NOT BE COPIED!!'    


                
