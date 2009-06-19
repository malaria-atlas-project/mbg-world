# example command line:
# run ECDISTRIBUTE_combineDistribExtractions r-095f3f60 extract_params_AF.py

# import libraries
import sys
from map_utils import amazon_ec
import numpy as np
from map_utils import checkAndBuildPaths

# deal with system arguments (expects two)
RESERVATIONID = str(sys.argv[1])  ## defines ID of reservation that contains the instances we will use on EC2
PARAMFILE = sys.argv[2]           ## defines name (inlduing '.py' extension) of parameter file to use

# initialise amazon S3 key object 
#S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

# set job distribution parameters
NINSTANCES = 1
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 1 # maximum number of tries before we give up on any individual job
STDOUTPATH = '/home/pwg/mbg-world/stdout_extraction/CombinedOutputSTDOUTERR_'+str(PARAMFILE.partition('.')[0])+'_'+str(time.ctime())+'/'
checkAndBuildPaths(STDOUTPATH,VERBOSE=True,BUILD=True)

# define files to upload to instance before any execution
UPLOADFILES=['/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']

# define any initialisation commands to exctue on instance before main job
INITCMDS=['bash /root/cloud_setup.sh']

# construct commands list
CMDS = ['"cd mbg-world/mbgw-scripts/;python extract_defineParameterFile.py '+str(PARAMFILE)+';python ECRUNSCRIPT_combineDistribExtractions.py True True True"'] 

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
