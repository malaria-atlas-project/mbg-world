## example call:
# run ECDISTRIBUTE_combineDistribExtractions r-49d0ae20 extract_params_AF.py
# run ECDISTRIBUTE_combineDistribExtractions r-ff8af396 extract_params_KE_eight.py
# run ECDISTRIBUTE_combineDistribExtractions r-cf8af3a6 extract_params_KE_nine.py
# run ECDISTRIBUTE_combineDistribExtractions r-75f28b1c extract_params_KE_ten.py
# run ECDISTRIBUTE_combineDistribExtractions r-75f28b1c extract_params_KE_eleven.py
# run ECDISTRIBUTE_combineDistribExtractions r-ff1c6796 extract_params_S1_twelve.py
# run ECDISTRIBUTE_combineDistribExtractions r-7d0b7014 extract_params_S2_thirteen.py
# run ECDISTRIBUTE_combineDistribExtractions r-7d0b7014 extract_params_S1_fourteen.py
# run ECDISTRIBUTE_combineDistribExtractions r-d1b6ccb8 extract_params_S1_fifteen.py
# run ECDISTRIBUTE_combineDistribExtractions r-85e2a6ec extract_params_S1_seventeen.py
# run ECDISTRIBUTE_combineDistribExtractions r-556a293c extract_params_AF.py
# run ECDISTRIBUTE_combineDistribExtractions r-5b6a2932 extract_params_AS1.py
# run ECDISTRIBUTE_combineDistribExtractions r-7bfbb712 extract_params_AM.py

# deal with system arguments (expects two)
## defines ID of reservation that contains the instances we will use on EC2
import sys
RESERVATIONID = str(sys.argv[1])
PARAMFILE = sys.argv[2]

# import libraries
from map_utils import amazon_ec
from map_utils import S3
import numpy as np
from map_utils import checkAndBuildPaths
import time

# initialise amazon S3 key object 
#S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

# set job distribution parameters
NINSTANCES = 1
MAXJOBSPERINSTANCE = 1
MAXJOBTRIES = 1 # maximum number of tries before we give up on any individual job
STDOUTPATH = '/home/pwg/mbg-world/stdout_extraction/CombinedOutputSTDOUTERR_'+str(PARAMFILE.partition('.')[0])+'_'+str(time.time())+'/'
checkAndBuildPaths(STDOUTPATH,VERBOSE=True,BUILD=True)

# define files to upload to instance before any execution
UPLOADFILES=['/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']

# define any initialisation commands to exctue on instance before main job
INITCMDS=['bash /root/cloud_setup.sh']

# construct commands list
CMDS = ['"cd mbg-world/mbgw-scripts/;python extract_defineParameterFile.py '+str(PARAMFILE)+';python ECRUNSCRIPT_combineDistribExtractions.py True True True"'] 

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
