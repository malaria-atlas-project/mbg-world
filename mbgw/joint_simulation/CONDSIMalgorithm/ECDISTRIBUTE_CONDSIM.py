# example command line:
# run ECDISTRIBUTE_CONDSIM r-933841fa CONDSIM_params_AF_eight.py 15
# run ECDISTRIBUTE_CONDSIM r-6f3b4206 CONDSIM_params_AF_nine.py 16
# run ECDISTRIBUTE_CONDSIM r-75f28b1c CONDSIM_params_AF_ten.py 16
# run ECDISTRIBUTE_CONDSIM r-6bcbb202 CONDSIM_params_AF_eleven.py 16
# run ECDISTRIBUTE_CONDSIM r-a9354ec0 CONDSIM_params_S1_twelve.py 16
# run ECDISTRIBUTE_CONDSIM r-9b2952f2 CONDSIM_params_S2_thirteen.py 16
# run ECDISTRIBUTE_CONDSIM r-c3176caa CONDSIM_params_S1_fourteen.py 16
# run ECDISTRIBUTE_CONDSIM r-dbd9a2b2 CONDSIM_params_S1_fifteen.py 16
# run ECDISTRIBUTE_CONDSIM r-77d3971e CONDSIM_params_S1_seventeen.py 16
# run ECDISTRIBUTE_CONDSIM r-81f4b0e8 CONDSIM_params_AF_eighteen.py 15
# python ECDISTRIBUTE_CONDSIM.py r-e1c68088 CONDSIM_params_AS1.py 15
# python ECDISTRIBUTE_CONDSIM.py r-87c583ee CONDSIM_params_AS2.py 15
# python ECDISTRIBUTE_CONDSIM.py r-df84c4b6 CONDSIM_params_AM.py 15

# import libraries
from map_utils import amazon_ec
from map_utils import S3
from map_utils import checkAndBuildPaths
import numpy as np
import time
import sys

#RESERVATIONID='r-e1eca288'
#PARAMFILE_PY='CONDSIM_params_AM.py'
#PARAMFILE_R=17
#NINSTANCES=1
#STAGE = 'BOTH'

# deal with system arguments (expects two)
RESERVATIONID = sys.argv[1]  ## defines ID of reservation that contains the instances we will use on EC2
PARAMFILE_PY = sys.argv[2]  ## defines name of python file housing the parmeter definitions (e.g. extract_params_AF.py)
PARAMFILE_R = int(sys.argv[3])  ## defines name of R file housing additoinal parmeter definitions for conditoinal simulation R scripts
NINSTANCES = int(sys.argv[4])
STAGE = sys.argv[5]

# initialise amazon S3 key object 
S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

STDOUTPATH = '/home/pwg/mbg-world/stdout_CONDSIM/DistributedOutputSTDOUTERR_'+str(PARAMFILE_PY.partition('.')[0])+'_'+str(time.ctime())+'/'
checkAndBuildPaths(STDOUTPATH,VERBOSE=True,BUILD=True)

if ((STAGE=='SETUP') | (STAGE == 'ALL')):

    # upload files
    print '\n*******************************'
    print 'STARTING UPLOADING FILES TO INSTANCES..'
    print '*******************************\n'
    MAXJOBSPERINSTANCE = 1
    MAXJOBTRIES = 2 #maximum number of tries before we give up on any individual job
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
    MAXJOBSPERINSTANCE = 1
    MAXJOBTRIES = 2 #maximum number of tries before we give up on any individual job
    UPLOADFILES=[]
    INITCMDS=[]
    CMDS = ['bash /root/cloud_setup.sh',]*NINSTANCES
    returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)
    CMDS = ['"cd /root/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm/;python CONDSIM_defineParameterFile.py '+str(PARAMFILE_PY)+';python ECRUNSCRIPT_CONDSIM_PREDOWNLOAD.py; python test_instance_setup.py"',]*NINSTANCES
    returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)  
    print '\n*******************************'
    print 'FINISHED EXECUTING INITILISATION COMMANDS ON INSTANCES..'
    print '*******************************\n'


if ((STAGE=='MAIN') | (STAGE == 'ALL')):

    # main jobs
    print '\n*******************************'
    print 'STARTING MAIN JOBS ON INSTANCES..'
    print '*******************************\n'
    MAXJOBSPERINSTANCE = 1
    MAXJOBTRIES = 1 #maximum number of tries before we give up on any individual job
    UPLOADFILES=[]
    INITCMDS=[]

    ## set realization number parameters
    n_total = 518#100 #600
    iter_per_job = 1
    NJOBS = n_total / iter_per_job
    STARTJOB = 378
    STOPJOB = 518 # this can be set to equal NJOBS, or a smaller number if we don;t want to do all NJOBS realisatios in one go - can continue with other realisations starting at i = STOPJOB

#    ##############TEMP
#    TEMPINDEX=np.array([951,950,931,925,923,917,914,908,891,890,87,879,877,869,868,858,848,838,836,82,818,803,780,76,765,762,752,74,743,721,706,69,692,690,65,659,641,57,53,4,467,465,436,42,426,41,415,3,399,397,38,37,360,358,34,340,333,320,31,316,314,311,308,299,295,283,27,279,277,275,266,264,260,230,228,214,212,205,202,18,17,173,170,16,165,157,156,14,13,137,136,11,119,107])
#    CMDS=[]
#    for i in xrange(STARTJOB,STOPJOB):
#        if (any(TEMPINDEX==i)): continue
#        CMDS.append('"cd /root/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm/;nice -n -20 python ECRUNSCRIPT_CONDSIM.py %i %i %i %i"'%(i,iter_per_job,NJOBS,PARAMFILE_R))
#    ##############TEMP


    ## construct main commands list
    CMDS = ['"cd /root/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm/;nice -n -20 python ECRUNSCRIPT_CONDSIM.py %i %i %i %i"'%(i,iter_per_job,NJOBS,PARAMFILE_R) for i in xrange(STARTJOB,STOPJOB)]

    ## finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
    startTime = time.time()
    returns = amazon_ec.map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=True,STDOUTPATH=STDOUTPATH)    
    endTime = time.time()-startTime
    print '\n*******************************'
    print 'FINISHED MAIN JOBS ON INSTANCES..'
    print '*******************************\n'
    print 'total run time for '+str(NJOBS)+' jobs with iter_per_job='+str(iter_per_job)+' on '+str(NINSTANCES)+' instances, with '+str(MAXJOBSPERINSTANCE)+' jobs per instance was: '+str(endTime)

