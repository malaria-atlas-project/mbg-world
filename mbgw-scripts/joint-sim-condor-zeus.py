import numpy as np
import pymc as pm
import mbgw
from mbgw.joint_simulation import *
import tables as tb
import mbgw.master_grid as mg
import os,sys
from boto_PYlib import *

print 'Imports done'

# Establish blocks based on spatial distance only.
i = int(sys.argv[1])
iter_per_job = int(sys.argv[2])
n_jobs = int(sys.argv[3])
region = sys.argv[4]
fname = sys.argv[5]
burn = int(sys.argv[6])
memmax = float(sys.argv[7])
thinning = int(sys.argv[8])
grid_lims = getattr(mg, region + '_lims')
nmonths = 288 
start_year = 1985
mask_name = 'st_mask5km-e_y-x+'
relp=1e-3

hf = tb.openFile(fname)
n_total = len(hf.root.chain0.PyMCsamples)
indices = np.array(np.linspace(burn, n_total, n_jobs+1), dtype=int)
my_start = indices[i]
my_end = indices[i+1]


ofdir = '/share/scratch/malaria-atlas-project/MAP-outputs/'
infile_base = fname.split('/')[-1].replace('.hdf5','')
outfile_base = 'realizations_mem_%i_%s.hdf5'%(memmax,'_'.join([infile_base, 'iterations', str(i*iter_per_job), str((i+1)*iter_per_job)]))
outfile_name = ofdir+outfile_base
 
print 'i: %i'%i
print 'iter_per_job: %i'%iter_per_job
print 'n_jobs: %i'%n_jobs
print 'region: %s'%region
print 'fname: %s'%fname
print 'burn: %i'%burn
print 'nmonths: %i'%nmonths
print 'start_year: %i'%start_year
print 'mask_name: %s'%mask_name
print 'relp: %f'%relp
print 'grid_lims: %s'%str(grid_lims)
print 'memmax: %i'%memmax
print 'Thinning: %i'%thinning

print 'Creating realizations'
create_many_realizations(my_start, iter_per_job, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, outfile_name, memmax, relp, mask_name, n_in_trace = my_end, thinning=thinning)
print 'Done'

print 'Uploading to boto'
S=S3('/home/oxg028/mbg-world/datafiles/s3code.txt')
S.uploadFileToBucket(infile_base.lower(),outfile_name,True,True)
print 'Done'
print 'Removing hdf5 archive locally'
os.remove(outfile_name)
print 'Done'