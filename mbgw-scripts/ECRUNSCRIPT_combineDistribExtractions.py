print "STARTING: ECRUNSCRIPT_combineDistribExtractions..\n"

# import libraries
from map_utils import checkAndBuildPaths
from extract_combineExtractions import *
from boto_PYlib import *
from extract_params import *

S=S3() # initialise key object

# deal with system arguments
BURDEN = True
PERPIXEL = True
PERCOUNTRY = True

if sys.argv[1] == 'False' : BURDEN=False
if sys.argv[2] == 'False' : PERPIXEL=False
if sys.argv[3] == 'False' : PERCOUNTRY=False

#print sys.argv[2]
#print type(sys.argv[2])
#print PERPIXEL
#print type(PERPIXEL)

if (PERPIXEL==True):

    # download from S3 contents of bucket 'distributedoutput_perpixel', will automatically build the local directory if necessary
    print '\n\tDownloading contents of S3 bucket "distributedoutput_perpixel" to local directory '+exportPathDistributed_perpixel
    S.downloadBucketContents('distributedoutput_perpixel',exportPathDistributed_perpixel,overwriteContent=False,VERBOSE=True)

    # build path for output to house combined per-pixel output maps
    print '\n\tChecking path for '+exportPathCombined_perpixel
    checkAndBuildPaths(exportPathCombined_perpixel,VERBOSE=True,BUILD=True)

    # download from S3 the other necessary files (optionally need 5km grump for burden map)
    print '\n\tDownloading lim5kmbnry file from S3..'
    S3bucketname = lim5kmbnry_path.split('/')[-2]
    print '\t\tS3bucketname: '+str(S3bucketname)
    S3filename = lim5kmbnry_path.split('/')[-1]
    print '\t\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,lim5kmbnry_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
    checkAndBuildPaths(lim5kmbnry_path,VERBOSE=True,BUILD=False)

    if (BURDEN==True):
        print '\n\tDownloading grump5km file from S3..'
        S3bucketname = grump5km_path.split('/')[-2]
        print '\t\tS3bucketname: '+str(S3bucketname)
        S3filename = grump5km_path.split('/')[-1]
        print '\t\tS3filename: '+str(S3filename)
        S.downloadFileFromBucket(S3bucketname,S3filename,grump5km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
        checkAndBuildPaths(grump5km_path,VERBOSE=True,BUILD=False)

    # now call extractSummaries_perpixel substituting in the formatted sys args 
    print '\n\tCalling combineDistribExtractions_perpixel'
    combineDistribExtractions_perpixel()

    # now upload the output back to the S3 storage

    print '\n\tuploading contents of exportPathCombined_perpixel to S3 bucket combinedoutput_perpixel'
    failCount = 0
    while failCount<=3:
        try:
            S.uploadDirectoryAsBucket('combinedoutput_perpixel',exportPathCombined_perpixel,uploadConstituentFiles=True,overwriteContent=True)
            break
        except RuntimeError:
            failCount+=1 
            if failCount<=3:
                print '\t\tuploading contents of exportPathCombined_perpixel to S3 bucket CombinedOutput_perpixel failed '+str(failCount) +' times: retrying..'
            else:
                print '\t\tuploading contents of exportPathCombined_perpixel to S3 bucket CombinedOutput_perpixel failed '+str(failCount) +' times: GIVING UP - FILE CONTENTS MAY NOT ALL HAVE COPIED!!'

if (PERCOUNTRY==True):

    # download from S3 contents of bucket 'distributedoutput_country', will automatically build the local directory if necessary
    print '\n\tDownloading contents of S3 bucket "distributedoutput_country" to local directory '+exportPathDistributed_perpixel
    S.downloadBucketContents('distributedoutput_country',exportPathDistributed_country,overwriteContent=False,VERBOSE=True)

    # download from S3 the salblim1km and salb1km file (used for calculating uniqueSalb and Nsalb / NsalbWholeCountries etc with function examineSalb)
    print '\n\tDownloading salblim1km file from S3..'
    S3bucketname = salblim1km_path.split('/')[-2]
    print '\t\tS3bucketname: '+str(S3bucketname)
    S3filename = salblim1km_path.split('/')[-1]
    print '\t\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,salblim1km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
    checkAndBuildPaths(salblim1km_path,VERBOSE=True,BUILD=False)

    print '\n\tDownloading salb1km file from S3..'
    S3bucketname = salb1km_path.split('/')[-2]
    print '\t\tS3bucketname: '+str(S3bucketname)
    S3filename = salb1km_path.split('/')[-1]
    print '\t\tS3filename: '+str(S3filename)
    S.downloadFileFromBucket(S3bucketname,S3filename,salb1km_path,overwriteContent=False,makeDirectory=True,VERBOSE=True)
    checkAndBuildPaths(salb1km_path,VERBOSE=True,BUILD=False)
    
    # build paths to directory to house uniqueSalb.txt, pixelN.txt,uniqueSalbAllCountries.txt, pixelNAllCountries.txt 
    checkAndBuildPaths(uniqueSalb_path.rpartition('/')[0],VERBOSE=True,BUILD=True)
    
    # run examineSalb on salblim1km to generate uniqueSalb.txt and pixelN.txt
    print '\n\trunning examineSalb on salblim1km..'
    temp=examineSalb (salblim1km_path,uniqueSalb_path,pixelN_path,ignore=np.array([-9999]))
    checkAndBuildPaths(uniqueSalb_path,VERBOSE=True,BUILD=False)
    checkAndBuildPaths(pixelN_path,VERBOSE=True,BUILD=False)

    # run examineSalb agauibn but this time using salb1km to generate uniqueSalbwholecountries.txt and pixelNwholecountries.txt - allows meanPR to be recalculated with all-country area as denominotor rather than just stable area
    print '\n\trunning examineSalb on salb1km..'
    temp=examineSalb (salb1km_path,uniqueSalbwholecountries_path,pixelNwholecountries_path,ignore=np.array([-9999]))
    checkAndBuildPaths(uniqueSalbwholecountries_path,VERBOSE=True,BUILD=False)
    checkAndBuildPaths(pixelNwholecountries_path,VERBOSE=True,BUILD=False)

    # build path for output to house combined per-pixel output maps
    print '\n\tChecking path for '+exportPathCombined_country
    checkAndBuildPaths(exportPathCombined_country,VERBOSE=True,BUILD=True)

    # now call extractSummaries_country substituting in the formatted sys args 
    print '\n\tCalling combineDistribExtractions_country'
    combineDistribExtractions_country()
            
    # now upload the main output back to the S3 storage            

    print '\n\tuploading contents of '+exportPathCombined_country+'to S3 bucket '+str(exportPathCombined_country.split('/')[-2].lower())
    failCount = 0
    while failCount<=3:
        try:
            S.uploadDirectoryAsBucket(exportPathCombined_country.split('/')[-2].lower(),exportPathCombined_country,uploadConstituentFiles=True,overwriteContent=True)
            break
        except RuntimeError:
            failCount+=1 
            if failCount<=3:
                print '\t\tuploading contents of '+exportPathCombined_country+' to S3 bucket '+exportPathCombined_country.split('/')[-2].lower()+' failed '+str(failCount) +' times: retrying..'
            else:
                print '\t\tuploading contents of '+exportPathCombined_country+' to S3 bucket '+exportPathCombined_country.split('/')[-2].lower()+' failed '+str(failCount) +' times: GIVING UP - FILE CONTENTS MAY NOT ALL HAVE COPIED!!'

    # now upload the examineSalb output back to the S3 storage            

    print '\n\tuploading contents of '+uniqueSalb_path.rpartition('/')[0]+'to S3 bucket misc'
    failCount = 0
    while failCount<=3:
        try:
            S.uploadDirectoryAsBucket('misc',uniqueSalb_path.rpartition('/')[0],uploadConstituentFiles=True,overwriteContent=True)
            break
        except RuntimeError:
            failCount+=1 
            if failCount<=3:
                print '\t\tuploading contents of '+uniqueSalb_path.rpartition('/')[0]+' to S3 bucket misc failed '+str(failCount) +' times: retrying..'
            else:
                print '\t\tuploading contents of '+uniqueSalb_path.rpartition('/')[0]+' to S3 bucket misc failed '+str(failCount) +' times: GIVING UP - FILE CONTENTS MAY NOT ALL HAVE COPIED!!'

print "FINISHED: ECRUNSCRIPT_combineDistribExtractions\n"

