# example command line:
# run ECDISTRIBUTE_CONDSIM r-933841fa CONDSIM_params_AF_eight.py 15
# run ECDISTRIBUTE_CONDSIM r-6f3b4206 CONDSIM_params_AF_nine.py 16
# run ECDISTRIBUTE_CONDSIM r-75f28b1c CONDSIM_params_AF_ten.py 16
# run ECDISTRIBUTE_CONDSIM r-6bcbb202 CONDSIM_params_AF_eleven.py 16
# run ECDISTRIBUTE_CONDSIM r-09c8b060 CONDSIM_params_S1_twelve.py 16


# import libraries
from map_utils import amazon_ec
from map_utils import S3
from map_utils import checkAndBuildPaths
import numpy as np
import time
import sys

# deal with system arguments (expects two)
RESERVATIONID = sys.argv[1]  ## defines ID of reservation that contains the instances we will use on EC2
PARAMFILE_PY = sys.argv[2]  ## defines name of python file housing the parmeter definitions (e.g. extract_params_AF.py)
PARAMFILE_R = int(sys.argv[3])  ## defines name of R file housing additoinal parmeter definitions for conditoinal simulation R scripts
#MAXJOBSPERINSTANCE = int(sys.argv[4]) 

# initialise amazon S3 key object 
S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

STDOUTPATH = '/home/pwg/mbg-world/stdout_CONDSIM/DistributedOutputSTDOUTERR_'+str(PARAMFILE_PY.partition('.')[0])+'_'+str(time.ctime())+'/'
checkAndBuildPaths(STDOUTPATH,VERBOSE=True,BUILD=True)

## upload files
NINSTANCES = 10
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 2 #maximum number of tries before we give up on any individual job
UPLOADFILES=['/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']
INITCMDS=[]
CMDS=[]
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    

## initialisation commands
NINSTANCES = 10
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 2 #maximum number of tries before we give up on any individual job
UPLOADFILES=[]
INITCMDS=[]
CMDS = ['bash /root/cloud_setup.sh','"cd /root/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm/;python CONDSIM_defineParameterFile.py '+str(PARAMFILE_PY)+';python ECRUNSCRIPT_CONDSIM_PREDOWNLOAD.py"',]*NINSTANCES
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    


## main jobs
NINSTANCES = 10
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job

# set realization number parameters
n_total = 80#100 #600
iter_per_job = 1
NJOBS = n_total / iter_per_job

STARTJOB = 0
STOPJOB = NJOBS # this can be set to equal NJOBS, or a smaller number if we don;t want to do all NJOBS realisatios in one go - can continue with other realisations starting at i = STOPJOB

# define files to upload to instance (from local machine) before any execution
UPLOADFILES=[]
INITCMDS=[]

# construct main commands list
CMDS = ['"cd /root/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm/;nice -n -20 python ECRUNSCRIPT_CONDSIM.py %i %i %i %i"'%(i,iter_per_job,NJOBS,PARAMFILE_R) for i in xrange(STARTJOB,STOPJOB)]

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
startTime = time.time()
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
endTime = time.time()-startTime

print 'total run time for '+str(NJOBS)+' jobs with iter_per_job='+str(iter_per_job)+' on '+str(NINSTANCES)+' instances, with '+str(MAXJOBSPERINSTANCE)+' jobs per instance was: '+str(endTime)

