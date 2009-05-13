# import libraries
from amazon_ec import *
from boto_PYlib import *
import numpy as np

# initialise amazon S3 key object 
S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

# define ID of reservation that contains the instances we will use on EC2
RESERVATIONID = 'r-5fe77536'

# set job distribution parameters
NINSTANCES = 1
MAXJOBSPERINSTANCE = 9
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
STDOUTPATH = '/home/pwg/mbg-world/extraction/DistributedOutput_perpixel/'

# set path to realisations on S3 and extract bucket and generic file name
realisations_path = '/mnt/qrypfpr010708_africa_run_9.10.2008_trial_two/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_FILESTARTREL_FILEENDREL.hdf5'
relBucket = realisations_path.rsplit('/')[-2]
relPath = realisations_path.rsplit('/')[-1]

# call queryRealizationsInBucket to obtain number and start/end realisation numbers of these realisation files
relDict = S.queryRealizationsInBucket(relBucket,relPath,VERBOSE=True)


# set realization number parameters
NRELS = relDict['Nrealisations']
NJOBS = relDict['Nfiles']

####################################TEMP
NJOBS = 35
####################################TEMP

FileStartRels = relDict['StartRelList']
FileEndRels = relDict['EndRelList']
NPER  = 2
NTOTALREL = NRELS*NPER

####################################TEMP
#NTOTALREL = 16
####################################TEMP

# construct main commands list
#CMDS = ['"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py %i %i %i %i None None True"'%(NPER,int(FileStartRels[i]),int(FileEndRels[i]),NTOTALREL) for i in xrange(NJOBS)]
CMDS=[]
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 100 101 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 14 15 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 22 23 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 28 29 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 32 33 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 37 38 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 44 45 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 48 49 18 None None True"')
CMDS.append('"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel.py 2 52 53 18 None None True"')

# define files to upload to instance (from local machine) before any execution
UPLOADFILES=['amazon_joint_sim.py','/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']

# define any initialisation commands to exctue on instance before main job
INITCMDS=['bash /root/cloud_setup.sh','"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perpixel_PREDOWNLOAD.py True"']

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
returns = map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
