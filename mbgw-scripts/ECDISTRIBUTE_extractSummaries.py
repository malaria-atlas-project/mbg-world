# example command line:
# run ECDISTRIBUTE_extractSummaries r-e57c028c extract_params_AF.py
# run ECDISTRIBUTE_extractSummaries r-ff8af396 extract_params_KE_eight.py
# run ECDISTRIBUTE_extractSummaries r-cf8af3a6 extract_params_KE_nine.py
# run ECDISTRIBUTE_extractSummaries r-75f28b1c extract_params_KE_ten.py
# run ECDISTRIBUTE_extractSummaries r-75f28b1c extract_params_KE_eleven.py
# run ECDISTRIBUTE_extractSummaries r-a9354ec0 extract_params_S1_twelve.py
# run ECDISTRIBUTE_extractSummaries r-9b2952f2 extract_params_S2_thirteen.py
# run ECDISTRIBUTE_extractSummaries r-c3176caa extract_params_S1_fourteen.py
# run ECDISTRIBUTE_extractSummaries r-d1b6ccb8 extract_params_S1_fifteen.py
# run ECDISTRIBUTE_extractSummaries r-77d3971e extract_params_S1_seventeen.py
# run ECDISTRIBUTE_extractSummaries r-f5f4b39c extract_params_AF.py
# python ECDISTRIBUTE_extractSummaries.py r-f5f4b39c extract_params_AS1.py
# python ECDISTRIBUTE_extractSummaries.py r-f5f4b39c extract_params_AS2.py
# python ECDISTRIBUTE_extractSummaries.py r-f5f4b39c extract_params_AM.py


# deal with system arguments (expects two)
import sys
RESERVATIONID = sys.argv[1]  ## defines ID of reservation that contains the instances we will use on EC2
PARAMFILE = sys.argv[2]  ## defines name of python file housing the parmeter definitions (e.g. extract_params_AF.py)

print 'importing local params from '+str(PARAMFILE.partition('.')[0])
localparams =__import__(PARAMFILE.partition('.')[0])

# import libraries
from map_utils import amazon_ec
from map_utils import S3
import numpy as np
from map_utils import checkAndBuildPaths
import time

# initialise amazon S3 key object 
S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

STDOUTPATH = '/home/pwg/mbg-world/stdout_extraction/DistributedOutputSTDOUTERR_'+str(PARAMFILE.partition('.')[0])+'_'+str(time.ctime())+'/'
checkAndBuildPaths(STDOUTPATH,VERBOSE=True,BUILD=True)

# set path to realisations on S3 and extract bucket and generic file name
relBucket = localparams.realisations_path.rsplit('/')[-2]
relPath = localparams.realisations_path.rsplit('/')[-1]

# call queryRealizationsInBucket to obtain number and start/end realisation numbers of these realisation files
relDict = S.queryRealizationsInBucket(relBucket,relPath,VERBOSE=True)

print '\nquerying bucket '+str(relBucket)+' : found '+str(relDict['Nrealisations'])+' realisations accross '+str(relDict['Nfiles'])+' files.'

# set realization number parameters
NRELS = relDict['Nrealisations']
NJOBS = relDict['Nfiles']

####################################TEMP
#NJOBS = 1
#NRELS = 1
####################################TEMP

FileStartRels = relDict['StartRelList']
FileEndRels = relDict['EndRelList']
NPER  = 10
NTOTALREL = NRELS*NPER

####################################TEMP 
#NTOTALREL = 1
####################################TEMP



# upload files
print '\n*******************************'
print 'STARTING UPLOADING FILES TO INSTANCES..'
print '*******************************\n'
NINSTANCES = 19
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
UPLOADFILES=['/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']
INITCMDS=[]
CMDS=[]
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
print '\n*******************************'
print 'FINISHED UPLOADING FILES TO INSTANCES..'
print '*******************************\n'


# initialisation commands
print '\n*******************************'
print 'STARTING EXECUTING INITILISATION COMMANDS ON INSTANCES..'
print '*******************************\n'
NINSTANCES = 19
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
UPLOADFILES=[]
INITCMDS=[]
CMDS = ['bash /root/cloud_setup.sh',]*NINSTANCES
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)
CMDS = ['"cd mbg-world/mbgw-scripts/;python extract_defineParameterFile.py '+str(PARAMFILE)+';python ECRUNSCRIPT_extractSummaries_PREDOWNLOAD.py True True True"',]*NINSTANCES
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)  
print '\n*******************************'
print 'FINISHED EXECUTING INITILISATION COMMANDS ON INSTANCES..'
print '*******************************\n'


# main jobs
print '\n*******************************'
print 'STARTING MAIN JOBS ON INSTANCES..'
print '*******************************\n'
NINSTANCES = 19
MAXJOBSPERINSTANCE = 2
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
UPLOADFILES=[]
INITCMDS=[]

## construct main commands list
CMDS = ['"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries.py %i %i %i %i None None True True True"'%(NPER,int(FileStartRels[i]),int(FileEndRels[i]),NTOTALREL) for i in xrange(NJOBS)]

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
startTime = time.time()
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=True,STDOUTPATH=STDOUTPATH)    
endTime = time.time()-startTime
print '\n*******************************'
print 'FINISHED MAIN JOBS ON INSTANCES..'
print '*******************************\n'
print 'total run time for '+str(NJOBS)+' jobs on '+str(NINSTANCES)+' instances, with '+str(MAXJOBSPERINSTANCE)+' jobs per instance was: '+str(endTime)







