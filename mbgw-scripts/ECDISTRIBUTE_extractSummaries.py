# example command line:
# run ECDISTRIBUTE_extractSummaries r-e57c028c extract_params_AF.py
# run ECDISTRIBUTE_extractSummaries r-ff8af396 extract_params_KE_eight.py
# run ECDISTRIBUTE_extractSummaries r-cf8af3a6 extract_params_KE_nine.py
# run ECDISTRIBUTE_extractSummaries r-75f28b1c extract_params_KE_ten.py
# run ECDISTRIBUTE_extractSummaries r-75f28b1c extract_params_KE_eleven.py

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

# set job distribution parameters
NINSTANCES = 1
MAXJOBSPERINSTANCE = 3
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
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

# define files to upload to instance (from local machine) before any execution
UPLOADFILES=['/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']

# define any initialisation commands to exctue on instance before main job
INITCMDS=['bash /root/cloud_setup.sh','"cd mbg-world/mbgw-scripts/;python extract_defineParameterFile.py '+str(PARAMFILE)+';python ECRUNSCRIPT_extractSummaries_PREDOWNLOAD.py True True True"']

####temp
#INITCMDS=['"cd mbg-world/mbgw-scripts/;python extract_defineParameterFile.py '+str(PARAMFILE)+';python ECRUNSCRIPT_extractSummaries_PREDOWNLOAD.py True True True"']
####temp


# construct main commands list
CMDS = ['"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries.py %i %i %i %i None None True True True"'%(NPER,int(FileStartRels[i]),int(FileEndRels[i]),NTOTALREL) for i in xrange(NJOBS)]

#from IPython.Debugger import Pdb
#Pdb(color_scheme='Linux').set_trace()   

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
startTime = time.time()
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
endTime = time.time()-startTime

print 'total run time for '+str(NJOBS)+' jobs with NPER='+str(NPER)+' on '+str(NINSTANCES)+' instances, with '+str(MAXJOBSPERINSTANCE)+' jobs per instance was: '+str(endTime)

