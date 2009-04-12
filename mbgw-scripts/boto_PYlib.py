# Author: Pete Gething
# Date: 5 March 2009
# License: Creative Commons BY-NC-SA
####################################

# set location of file containng keys to S3
#keyPath = '/root/s3code.txt'
from extract_params import keyPath

# import libraries
import boto
from boto.s3.connection import S3Connection
import random
import sys
import string
import os
import time
import md5
from socket import gethostname
import numpy as np
from map_utils import checkAndBuildPaths

##############################################################################################################################
def S3connectKeys(keyNo):

    a=file(keyPath).read()
    key = a.split(',')[keyNo]
    key=key.strip("\'")
    key=key.strip(" \'")
    return(key)
##############################################################################################################################
def getFileAgeMinutes(filekey):

    # get string date this file was last modified on S3
    lm = filekey.last_modified

    # convert this string into year,month,day,hr,mn,sec values to match time.gmtime
    monthList=list(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
    monthDay=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    lmlist=lm.split()
    lm_year = int(lmlist[3])
    lm_day = int(lmlist[1])
    lm_mon = monthList.index(lmlist[2])
    lm_yday = monthDay.cumsum()[lm_mon-1] + lm_day
    timestr1 = lmlist[4]
    timestr2 = timestr1.partition(':')
    timestr3 = timestr2[2].partition(':')
    lm_hour = int(timestr2[0])
    lm_min = int(timestr3[0])
    lm_sec = int(timestr3[2])
    timeModSec = (3600*24*365*lm_year) + (3600*24*lm_yday) + (3600*lm_hour) + (60*lm_min) + (lm_sec)

    # compare time now with time file modified to get file age
    timeNow = time.gmtime()
    timeNowSec = (3600*24*365*timeNow[0]) + (3600*24*timeNow[7]) + (3600*timeNow[3]) + (60*timeNow[4]) + (timeNow[5])
    ageSec = timeNowSec - timeModSec
    ageMin =ageSec/60.

    return(ageMin)
##############################################################################################################################
def uploadDirectoryAsBucket(bucketName,directoryPath,uploadConstituentFiles,overwriteContent):

    '''
    Allows simple copying of a local directory (and optionally its constituesnt files)to S3 as
    a bucket of the same name containing files of the same name.
    
    params to pass:
    bucketName              : (string) name to give bucket when it is created
    directoryPath           : (string) full path of directory on local machine that will be copied
    uploadConstituentFiles  : (logical) do we want to upload the files in this directory as objects in the bucket?
    overwriteContent        : (logical) if a file already exists in the bucket do we want to overwrite it? if False then only new files in this directory will be copied to bucket
    '''
    
    # run checks on this directory path
    if (os.path.exists(directoryPath) != True) :
        print "ERROR!!! path "+directoryPath+" does not exist: EXITING!!!\n"
        return(-9999)

    if (os.path.isdir(directoryPath) != True) :
        print "ERROR!!! path "+directoryPath+" is not a directory: EXITING!!!\n"
        return(-9999)

    # make sure path does not end with a '/' to ensure split works
    if directoryPath[-1] == '/':
        directoryPath = directoryPath[slice(0,len(directoryPath)-1,1)]

    # seperate directory name from path
    pathsplit = os.path.split(directoryPath)
    Path = pathsplit[0]
    directoryName = pathsplit[1]

    # establish connection to S3 account
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))

    # create the bucket
    bucket = conn.create_bucket(str(bucketName))

    # upload a metadata file containing the path of the origin file
    k = boto.s3.key.Key(bucket)
    k.key='bucketOrigin.txt'
    k.set_contents_from_string('hostname:'+str(gethostname())+' path:'+str(directoryPath))

    # optionally, upload all files in this directory to this bucket
    if uploadConstituentFiles==True:

        for fname in os.listdir(directoryPath):
        
            # define full local path+filename for this file
            filename = directoryPath+'/'+str(fname)

            # check this is actually a file and not a directory before continuing
            if (os.path.isfile(filename) != True) : continue
            
            # if we are not allowing overwriting of files in this bucket, check if this file already exists on S3 and if so skip overwrite
            if overwriteContent==False:
                if (CheckFileExistsInBucket(directoryName,fname) == True) : continue

            # establish key object
            k = boto.s3.key.Key(bucket)

            # give this key the same name as the file
            k.key = str(fname)

            # pass this file to S3 using the key
            k.set_contents_from_filename(filename)

            # finally, check that this file made it to S3
            md5string = md5.new(file(filename).read()).hexdigest()
            if (CheckFileExistsInBucket(bucketName,fname,md5check = md5string) != True) :
                print 'ERROR!! final check revealed file "'+str(filename)+'" did not copy succesfully to S3 bucket "'+str(bucketName)+'"'
##############################################################################################################################
def makeEmptyBucket(bucketName):

    '''
    Simply creates a new empty bucket with the specified name, checking that
    no bucket of the same name already exists.
    
    params to pass:
    bucketName              :(string) name to give bucket when it is created
    '''
 
    #ensure input is string
    bucketName=str(bucketName)

    # establish connection to S3 account and to specified bucket
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))    

    # check this bucket does not already exist
    if (conn.lookup(bucketName) is not None):
        print 'WARNING!!! requested bucket "'+str(bucketName)+'" already exists on S3 !!!'
        return(-9999)
    
    # create bucket with given name
    bucket = conn.create_bucket(bucketName)
    
    # check bucket now exists
    if (conn.lookup(bucketName) is None):
        print 'WARNING!!! requested bucket "'+str(bucketName)+'" does not appear to have been made on  S3 !!!'
        return(-9999)    
    
    return(0)
##############################################################################################################################
def uploadFileToBucket(bucketName,filePath,overwriteContent,makeBucket,VERBOSE=True):

    '''
    Allows simple copying of a single local file to a specified bucket.
    
    params to pass:
    bucketName              : (string) name of bucket to send file
    filePath                : (string) full path of file on local machine that will be copied
    overwriteContent        : (logical) if an object of this name already exists in the bucket, do we want to overwrite?
    makeBucket              : (logical) if the specified bucket does not exist, do we want to make it?
    '''

    #ensure input is string
    bucketName=str(bucketName)
    filePath=str(filePath)

    # check this is actually a file and it exists before continuing
    if (os.path.exists(filePath) != True) :
        print 'ERROR!!! requested file "'+str()+'" cannot be found: EXITING!!'
        return(-9999)

    if (os.path.isfile(filePath) != True) :
        print 'ERROR!!! requested file "'+str()+'" to upload to bucket is not a file: EXITING!!'
        return(-9999)
    
    # establish connection to S3 account 
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))        

    # check sepcified bucket  exists and optionally make it
    if (conn.lookup(bucketName) is None):
        if makeBucket==False:
            if VERBOSE==True: print 'WARNING!!! requested bucket "'+str(bucketName)+'" does not exist on  S3 and makeBucket==False!!!'
            return(-9999)

        if makeBucket==True:
            temp=makeEmptyBucket(bucketName)
            if(temp==-9999):
                print 'ERROR!!! recieved error return when trying to make new bucket "'+str(bucketName)+'" EXITING!!!'


    # obtain filename from path
    fileName = os.path.split(filePath)[1]
        
    # if worried about overwriting check if file exists already and abort if it does
    if overwriteContent==False:
        if (CheckFileExistsInBucket(bucketName,fileName) == True) :
            print 'WARNING!!! file "'+str(fileName)+'" already exists in bucket "'+str(bucketName)+'" and overwriteContent==False: EXITING!!'
            return(-9999)

    # now go ahead and upload file to the bucket

    ## link to the bucket with same name as directory
    bucket = conn.get_bucket(bucketName)

    ## establish key object
    k = boto.s3.key.Key(bucket)

    ## give this key the same name as the file
    k.key = str(fileName)

    ## pass this file to S3 using the key
    k.set_contents_from_filename(filePath)
    
    # finally, check that this file made it to S3
    md5string = md5.new(file(filePath).read()).hexdigest()
    if (CheckFileExistsInBucket(bucketName,fileName,md5check = md5string) != True) :
        print 'ERROR!! final check revealed file "'+str(filePath)+'" did not copy succesfully to S3 bucket '+str(bucketName)

    return(0)
##############################################################################################################################
def isLOCALFILEIdenticalToS3FILE(bucketName,fileNameInBucket,localFilePath):

    '''
    Checks whther a local file and a file on S3 are identical, according to their md5 strings. Does all the necessary checks
    and returns a True/False accordingly
    
    params to pass:
    bucketName              : (string) name of bucket that file of interest is located in 
    fileNameInBucket        : (string) name of file of interest in the bucket
    localFilePath           : (logical) full path to local file of interest
    '''

    # check local file exists
    if(checkAndBuildPaths(localFilePath,VERBOSE=True,BUILD=False)==-9999): return(False)
    
    # get md5 string for local file
    md5string = md5.new(file(localFilePath).read()).hexdigest()
    
    # get md5 sring for file in bucket    
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))

    ## check this bucket exits
    if (conn.lookup(bucketName) is None):
        print 'WARNING!!! requested bucket "'+str(bucketName)+'" does not exist on S3 !!!'
        return(False)

    ## check the file exists on this bucket    
    if(CheckFileExistsInBucket(bucketName,fileNameInBucket,VERBOSE=True)!=True): return(False)  

    ## get md5 string for this file in the bucket
    bucket = conn.get_bucket(bucketName)
    filekey=bucket.get_key(fileNameInBucket)
    md5_s3 = filekey.etag.strip(str('"'))

    # compare this to the passed md5 string
    if (md5string != md5_s3):
        print 'ERROR!!! the md5 string of file "'+str(fileNameInBucket)+'" in bucket "'+str(bucketName)+'" does not match that of local file at "'+str(localFilePath)+'"  !!!!'
        return(False)

    # if these tests passed, then return True
    return(True)    
##############################################################################################################################
def CheckFileExistsInBucket(bucketName,fileNameInBucket, md5check=None, maxAgeMins=None, VERBOSE=False):

    '''
    Checks whether a file of the specified name exists in the specified bucket. Various options:
    1. md5check=None, maxAgeMins=None  = simply checks by filename
    2. maxAgeMins=float(minutes)  = checks that file exists and is younger than maxAgeMins
    2. md5check=str(md5str)  = checks that file exists and has identical md5 string to that passed by md5check
    
    params to pass:
    bucketName              : (string) name of bucket that file of interest is located in 
    fileNameInBucket        : (string) name of file of interest in the bucket
    md5check                : (string) an md5 string - if passed, will check that file in bucket has identical string
    maxAgeMins                : (float) if passed, will check that file in bucket is younger (since last modified) than this age in minutes
    '''
    
    # establish connection to S3 account and to specified bucket
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))
    
    # check this bucket exits
    if (conn.lookup(bucketName) is None):
        print 'WARNING!!! requested bucket "'+str(bucketName)+'" does not exist on S3 OR the connection failed using these access keys!!!'
        return(False)
    
    # check this file exists within this bucket
    bucket = conn.get_bucket(bucketName)
    keylist=[]
    rs=bucket.list()
    for key in rs:
        keylist.append(str(key.name))
    
    if (keylist.count(fileNameInBucket)==0) :
        if VERBOSE==True: print 'WARNING!!! file "'+str(fileNameInBucket)+'" not found in bucket "'+str(bucketName)+'"'
        return(False)
        
    # optionally check that the file was updated within the specified time
    if maxAgeMins is not None:
    
        filekey=bucket.get_key(fileNameInBucket)
        fileAge = getFileAgeMinutes(filekey)
        
        if fileAge>maxAgeMins:
            if VERBOSE==True: print 'WARNING!!! file "'+str(fileNameInBucket)+'" is older ('+str(fileAge)+' mins) than maxAgeMins ('+str(maxAgeMins)+' mins)'
            return(False)

    # optionally check that the file in the bucket has the same md5 string as that passed (e.g originating from a local file to check is identical)
    if md5check is not None:
       
        # get md5 string for this file in the bucket
        filekey=bucket.get_key(fileNameInBucket)
        md5_s3 = filekey.etag.strip(str('"'))

        # copmare this to the passed md5 string
        if (md5check != md5_s3):
            print 'ERROR!!! the md5 string of file "'+str(fileNameInBucket)+'" in bucket "'+str(bucketName)+'" does not match that passed !!!!'
            return(False)
    
    # if these tests passed, then return True
    return(True)
##############################################################################################################################
def downloadFileFromBucket(bucketName,fileNameInBucket,filePathAtDestination,overwriteContent,makeDirectory,VERBOSE=True):

    '''
    Downloads a single file from a bucket to a local directory
    
    params to pass:
    bucketName              : (string) name of bucket in which file of interest is lcoated
    fileNameInBucket        : (string) name of file of interest in bucket
    filePathAtDestination   : (string) full path - including filename itself - of the file once it has been downloaded
    overwriteContent        : (logical) if the specified target file already exists locally, do we overwrite?
    makeDirectory           : (logical) if the specified target file path includes new directories, should we make these?
    '''

    # check file exists in bucket
    if (CheckFileExistsInBucket(bucketName,fileNameInBucket,VERBOSE=True) != True) : return(-9999)
        
    # check destination directory exists and optionally build if not
    if(checkAndBuildPaths (os.path.split(filePathAtDestination)[0],VERBOSE=False,BUILD=makeDirectory)==-9999): return(-9999)

    # if we are not overwriting, then check if file ealrady exists locally and abort if it does
    if overwriteContent==False:
        if (os.path.exists(filePathAtDestination)==True):
            print 'ERROR!! file "'+str(filePathAtDestination)+'" already exists and overwriteContent==False: Exiting !!!' 
            return(-9999)
            
    # if we are overwriting, check that existing file is not actually a directory
    if overwriteContent==True:
        if(os.path.isdir(filePathAtDestination)==True):
            print 'ERROR!! a directory ("'+str(fileNameInBucket)+'") exists at path "'+str(filePathAtDestination)+'" with same name as file trying to download: EXITING!!'
            return(-9999)

    # establish connection to S3 account and to specified bucket
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))
    bucket = conn.get_bucket(bucketName)

    # establish key object
    filekey=bucket.get_key(fileNameInBucket)

    # pass the contents of file on s3 to the local file
    filekey.get_contents_to_filename(filePathAtDestination)

    # finally, check that this file made it from S3 to local destination

    ## first check there is even a file of this name at local destination
    if(os.path.exists(filePathAtDestination)!= True):
        print 'ERROR!! final check revealed file "'+str(filePathAtDestination)+'" did not copy succesfully from S3 file "'+str(fileNameInBucket)+'" in bucket "'+str(bucketName)+'"'
        return(-9999)

    ## then check the md5 keys match
    md5_s3 = filekey.etag.strip(str('"'))
    md5string = md5.new(file(filePathAtDestination).read()).hexdigest()

    if(md5string != md5_s3):
        print 'ERROR!! final check revealed file "'+str(filePathAtDestination)+'" did not copy succesfully from S3 file "'+str(fileNameInBucket)+'" in bucket "'+str(bucketName)+'"'
        return(-9999)
##############################################################################################################################
def downloadBucketContents(bucketName,targetDirectoryPath,overwriteContent,VERBOSE=True):

    '''
    Downloads contents of an entire bucket to a specified local file.
    
    params to pass:
    bucketName              : (string) name of bucket of interest  
    targetDirectoryPath     : (string) path to target directory. If this includes new directories, these will be built if possible
    overwriteContent        : (logical) if the specifid pathis to an existing directory, and there are existing files with the same name as those in the bucket, do we overwrite?
    '''
    
    # check bucket exists on S3
    conn = S3Connection(S3connectKeys(0), S3connectKeys(1))
    if (conn.lookup(bucketName) is None):
        print 'WARNING!!! requested bucket "'+str(bucketName)+'" does not exist on S3 !!!'
        return(-9999)
        
    # check target local directory exists and if not then build it        
    if(checkAndBuildPaths(targetDirectoryPath,BUILD=True)==-9999):
        print 'ERROR!! problem building target directory "'+str(targetDirectoryPath)+'" : EXITING!!!'
        return(-9999)

    # get list of files already in target directory
    existinglocalfiles = os.listdir(targetDirectoryPath)

    # loop through all files in the bucket
    bucket = conn.get_bucket(bucketName)
    rs=bucket.list()
    for key in rs:
    
        # if not overwriting, check no file exists in local directory with same name as this file
        if overwriteContent==False:
            if (existinglocalfiles.count(str(key.name))>0):
                if VERBOSE==True: print 'WARNING!!! file "'+str(key.name)+'" already present in local directory "'+str(targetDirectoryPath)+'" and overwriteContent==False '
                continue

        # if we are overwriting, check that existing file is not actually a directory
        if overwriteContent==True:
            if(os.path.isdir(targetDirectoryPath+str(key.name))==True):
                print 'ERROR!! a directory ("'+str(key.name)+'") exists at path "'+str(targetDirectoryPath)+'" with same name as file trying to download: EXITING!!'
                return(-9999)

        
        # build full target filepath
        if targetDirectoryPath[-1] != '/': targetDirectoryPath = targetDirectoryPath+'/'
        filePathAtDestination = targetDirectoryPath+str(key.name)

        # now copy this file from S3 bucket to local directory
        key.get_contents_to_filename(filePathAtDestination)        
        
        # check file has made it to destination
        
        ## first check there is even a file of this name at local destination
        if(os.path.exists(filePathAtDestination)!= True):
            print 'ERROR!! final check revealed file "'+str(filePathAtDestination)+'" did not copy succesfully from S3 file "'+str(key.name)+'" in bucket "'+str(bucketName)+'"'
            return(-9999)

        ## then check the md5 keys match
        md5_s3 = key.etag.strip(str('"'))
        md5string = md5.new(file(filePathAtDestination).read()).hexdigest()

        if(md5string != md5_s3):
            print 'ERROR!! final check revealed file "'+str(filePathAtDestination)+'" did not copy succesfully from S3 file "'+str(key.name)+'" in bucket "'+str(bucketName)+'"'
            return(-9999)
##############################################################################################################################