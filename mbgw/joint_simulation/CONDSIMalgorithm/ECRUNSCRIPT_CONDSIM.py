print 'Starting ECRUNSCRIPT_CONDSIM....'

import numpy as np
import pymc as pm
import mbgw
import time
from mbgw.joint_simulation import *
import tables as tb
import mbgw.master_grid as mg
import os,sys
from map_utils import S3
from CONDSIM_params import *

print 'Imports done'


# handle system parameters
i = int(sys.argv[1])
iter_per_job = int(sys.argv[2])
n_jobs = int(sys.argv[3])
paramfileINDEX = sys.argv[4]

# define some derived parameters

## file paths and grid set up
grid_lims = getattr(mg, region + '_lims')
mask_name = lim5kmbnry_path.split('/')[-1].replace('.hdf5','')
hf = tb.openFile(trace_path)
infile_base = trace_path.split('/')[-1].replace('.hdf5','')

#######################TEMP
outfile_base = 'realizations_mem_%i_%s.hdf5'%(memmax,'_'.join([infile_base, 'iterations', str(i*iter_per_job), str((i+1)*iter_per_job)]))
#outfile_base = 'realizations_mem_%i_%s.hdf5'%(memmax,'_'.join([infile_base, 'iterationsCOND', str(i*iter_per_job), str((i+1)*iter_per_job)]))
#######################TEMP
outfile_name = realizations_path+outfile_base

## iterations and indices set up
n_total = len(hf.root.chain0.PyMCsamples)

if (burn+n_jobs)>n_total:
    raise ValueError ('insufficient realisations in trace ('+str(n_total)+') to allow burn = '+str(burn)+' and n_jobs= '+str(n_jobs))

indices = np.array(np.linspace(burn, n_total, n_jobs+1), dtype=int)
my_start = indices[i]
my_end = indices[i+1]

 
print 'i: %i'%i
print 'iter_per_job: %i'%iter_per_job
print 'n_jobs: %i'%n_jobs
print 'region: %s'%region
print 'burn: %i'%burn
print 'nmonths: %i'%nmonths
print 'start_year: %i'%start_year
print 'mask_name: %s'%mask_name
print 'relp: %f'%relp
print 'grid_lims: %s'%str(grid_lims)
print 'memmax: %i'%memmax
print 'Thinning: %i'%thinning

print 'Creating realizations'

if region == "AM":
    print '\nrunning preliminary test square for AMS1 region'
    grid_lims_AMS1 = getattr(mg, "AMS1" + '_lims')
    create_many_realizations(my_start, iter_per_job, hf.root.chain0, hf.root.metadata, grid_lims_AMS1, start_year, nmonths, realizations_path+'TESTSQUARE'+outfile_base, memmax, relp, mask_name, n_in_trace = my_end, thinning=thinning,paramfileINDEX=paramfileINDEX,NinThinnedBlock=NinThinnedBlock,merged_urb=merged_urb,TESTRANGE=TESTRANGE,TESTSQUARE=True) 
    print '\nDONE preliminary test square for AMS1 region...will go on to full version'


t1=time.time()
create_many_realizations(my_start, iter_per_job, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, outfile_name, memmax, relp, mask_name, n_in_trace = my_end, thinning=thinning,paramfileINDEX=paramfileINDEX,NinThinnedBlock=NinThinnedBlock,merged_urb=merged_urb,TESTRANGE=TESTRANGE)

print 'Total time for realizations was '+str(time.time()-t1)

#from IPython.Debugger import Pdb
#Pdb(color_scheme='Linux').set_trace()

print 'Uploading output realization :'+str(outfile_name)+' to S3 bucket : '+str(realizations_path.split('/')[2])
S=S3(keyPath)
S.uploadFileToBucket(realizations_path.split('/')[2],outfile_name,True,True)

print 'Finished ECRUNSCRIPT_CONDSIM'







