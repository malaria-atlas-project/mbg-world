# example command line:
# run ECDISTRIBUTE_extractSummaries r-41422228 extract_params_AF.py

# deal with system arguments (expects two)
import sys
RESERVATIONID = sys.argv[1]  ## defines ID of reservation that contains the instances we will use on EC2
PARAMFILE = sys.argv[2]  ## defines name of python file housing the parmeter definitions (e.g. extract_params_AF.py)

print 'importing local params from '+str(PARAMFILE.partition('.')[0])
localparams =__import__(PARAMFILE.partition('.')[0])

# import libraries
from amazon_ec import * 
from boto_PYlib import *
import numpy as np
from map_utils import checkAndBuildPaths
import time

# initialise amazon S3 key object 
S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

# set job distribution parameters
NINSTANCES = 20
MAXJOBSPERINSTANCE = 3
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
STDOUTPATH = '/home/pwg/mbg-world/extraction/DistributedOutputSTDOUTERR_'+str(PARAMFILE.partition('.')[0])+'_'+str(time.ctime())+'/'
checkAndBuildPaths(STDOUTPATH,VERBOSE=True,BUILD=True)

# set path to realisations on S3 and extract bucket and generic file name
relBucket = localparams.realisations_path.rsplit('/')[-2]
relPath = localparams.realisations_path.rsplit('/')[-1]

# call queryRealizationsInBucket to obtain number and start/end realisation numbers of these realisation files
relDict = S.queryRealizationsInBucket(relBucket,relPath,VERBOSE=True)

# set realization number parameters
NRELS = relDict['Nrealisations']
NJOBS = relDict['Nfiles']

####################################TEMP
#NJOBS = 237
#NRELS = 180
####################################TEMP

FileStartRels = relDict['StartRelList']
FileEndRels = relDict['EndRelList']
NPER  = 100
NTOTALREL = NRELS*NPER

####################################TEMP 
#NTOTALREL = 2
####################################TEMP

# define files to upload to instance (from local machine) before any execution
UPLOADFILES=['amazon_joint_sim.py','/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']

# define any initialisation commands to exctue on instance before main job
INITCMDS=['bash /root/cloud_setup.sh','"cd mbg-world/mbgw-scripts/;python extract_defineParameterFile.py '+str(PARAMFILE)+';python ECRUNSCRIPT_extractSummaries_PREDOWNLOAD.py False True False"']

# construct main commands list
CMDS = ['"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries.py %i %i %i %i None None False True False"'%(NPER,int(FileStartRels[i]),int(FileEndRels[i]),NTOTALREL) for i in xrange(NJOBS)]

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
startTime = time.time()
returns = map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=True,STDOUTPATH=STDOUTPATH)    
endTime = time.time()-startTime

print 'total run time for '+str(NJOBS)+' jobs with NPER='+str(NPER)+' on '+str(NINSTANCES)+' instances, with '+str(MAXJOBSPERINSTANCE)+' jobs per instance was: '+str(endTime)

