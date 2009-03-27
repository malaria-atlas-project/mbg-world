import sys
import time
import subprocess
import numpy as np
import tables as tb

#####
a=time.time()
####


# sepcify some parameters
STARTREL = 0        # first full realisation to include
ENDREL = 8        # last full realisation to include (does up to but not including this index)
NPER = 50           # how many draws of the nugget, per full realisation
RELPERJOB = 1       # ideally, how many realisations do we want to ascribe to each job?
MAXJOBS = 8         # how many processes do we want running at a time
WAITINTERVAL = 30  # wait interval in seconds
VERBOSE = True      # messages or not
STDOUTPUT = '/home/pwg/mbg-world/extraction/DistributedOutput/scriptoutput/'   # if string given, this is where the standard output and standard error will be sent, if 0, then none sent

###########################################################################################################################
def returnNjobsRunning(splist):

    # two simple functions to allow 'x.Poll' and 'x is None' to be applied accross elements of a list
    def applyPoll(sp):
        return(sp.poll())

    def apply_is_None(sp):
        return(sp is None) 

    # apply the subprocess method .Poll to each element of the subprocess list to get list of statuses 
    jobstatus = np.array(map(applyPoll,splist))

    # query how many elements of this list are 'None' which Poll returns if job is still running
    sumRunning = sum(map(apply_is_None,jobstatus))

    # return this sum
    return({'sumRunning':sumRunning,'jobstatus':jobstatus})

###########################################################################################################################


# initialise list to hold subprocess objects (containing info about each subprocess)
splist = []

# initialise counter for output messages
if VERBOSE: cnt = 0

# import link to hdf5 realisation block in order to query how many realisations are present - to check if these tally with those asked for above..
from extract_params import filename
from map_utils import checkAndBuildPaths
checkAndBuildPaths(filename,VERBOSE=False,BUILD=False)
hf = tb.openFile(filename) 
RelsInBlock = hf.root.realizations.shape[0]

# run check that number of realisations asked for is actually present in block
NREL = ENDREL - STARTREL
if(NREL>RelsInBlock):
    print "ERROR!! have asked for "+str(NREL)+" realisations but only "+str(RelsInBlock)+ "present in hdf5 block: EXITING!!"
    raise ValueError

# define LUT of start and end realisation for each job
startRels = np.arange(STARTREL,ENDREL,RELPERJOB)
endRels = startRels + RELPERJOB
endRels[-1] = ENDREL
NJOBS = len(startRels)

# set off jobs
for ii in xrange(NJOBS):

    startRel = startRels[ii]
    endRel = endRels[ii]

    # if asked for, define stndout and stderr files for this job
    stdout=None
    stderr=None  
    if STDOUTPUT != 0:
        stdout = file(STDOUTPUT+'stdout_r'+str(startRel)+"to"+str(endRel-1)+".txt",'w')
        stderr = file(STDOUTPUT+'stderr_r'+str(startRel)+"to"+str(endRel-1)+".txt",'w')
        
     # define command line syntax for this job    
    args = (['python','/home/pwg/mbg-world/mbgw-scripts/extract_PYlib.py',str(NPER),str(startRel),str(endRel)])

    # start this subprocess, and assign the resulting subprocess object to the splist
    splist.append(subprocess.Popen(args=args,stdout=stdout,stderr=stderr))
    if VERBOSE:
        print 'sent job '+str(cnt)
        cnt=cnt+1

    # if number of jobs running has reached MAXJOBS, keep checking for jobs to finish before sending more
    while ((returnNjobsRunning(splist)['sumRunning'] >= MAXJOBS) & (returnNjobsRunning(splist)['sumRunning']<NJOBS)):
        if VERBOSE: print("waiting to start more jobs at time "+ str(time.time()))
        time.sleep(WAITINTERVAL)

# once all jobs have been dispatched, wait for all to finish before continuing
while returnNjobsRunning(splist)['sumRunning'] !=0:
    if VERBOSE: print("waiting for all jobs to finish at time "+ str(time.time()))
    time.sleep(WAITINTERVAL)

# check for any jobs that did not execute successfully
sumSuccess = sum(np.array(returnNjobsRunning(splist)['jobstatus'])==0)
sumFail = sum(np.array(returnNjobsRunning(splist)['jobstatus'])>0)
if VERBOSE: print 'successful jobs: '+str(sumSuccess)+' of '+str(len(splist))
if sumFail!=0: print 'failed jobs: '+str(sumFail)+' of '+str(len(splist))

# run script to recombine distruted output
from extract_combineDistribExtractions import *

#######
print("TOTAL TIME: "+(str(time.time()-a)))
######


