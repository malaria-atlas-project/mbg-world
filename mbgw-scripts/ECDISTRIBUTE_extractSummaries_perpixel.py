# import libraries
from amazon_ec import *
from boto_PYlib import *
import numpy as np

# initialise amazon S3 key object 
#S=S3(keyPath='/home/pwg/mbg-world/mbgw-scripts/s3code.txt')

# define ID of reservation that contains the instances we will use on EC2
RESERVATIONID = 'r-bc6ffed5'

# set job distribution parameters
NINSTANCES = 2
MAXJOBSPERINSTANCE = 2
MAXJOBTRIES = 3 #maximum number of tries before we give up on any individual job
STDOUTPATH = '/home/pwg/mbg-world/extraction/DistributedOutput_perpixel/'

# set realization number parameters
NRELS = 3
NPER  = 2
NRELSPERJOB = 1.
FIRSTREL = 0

NTOTALREL = NRELS*NPER
NJOBS = np.ceil(NRELS/NRELSPERJOB)

FileStartRels = FIRSTREL + np.arange(NJOBS)*NRELSPERJOB
FileEndRels = FileStartRels+NRELSPERJOB

# construct commands list
CMDS = ['"cd mbg-world/mbgw-scripts/;python ECRUNSCRIPT_extractSummaries_perPixel.py %i %i %i %i None None True"'%(NPER,FileStartRels[i],FileEndRels[i],NTOTALREL) for i in xrange(NJOBS)]

# define files to upload to instance before any execution
UPLOADFILES=['amazon_joint_sim.py','/home/pwg/mbg-world/mbgw-scripts/cloud_setup.sh','/home/pwg/mbg-world/mbgw-scripts/s3code.txt']

# define any initialisation commands to exctue on instance before main job
INITCMDS=['bash /root/cloud_setup.sh']

# finally, call local function map_jobs from amazon_ec module to distribute these jobs on EC2
returns = map_jobs(RESERVATIONID,NINSTANCES,MAXJOBSPERINSTANCE,MAXJOBTRIES,cmds=CMDS, init_cmds=INITCMDS,upload_files=UPLOADFILES, interval=20,shutdown=False,STDOUTPATH=STDOUTPATH)    
