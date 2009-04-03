import numpy as np
import pymc as pm
import mbgw
from mbgw.joint_simulation import *
import tables as tb
import mbgw.master_grid as mg
import sys

print 'Imports done'

# Establish blocks based on spatial distance only.
n_blocks_x = 10
n_blocks_y = 10
N_nearest = 300

i = int(sys.argv[1])
iter_per_job = int(sys.argv[2])
n_jobs = int(sys.argv[3])
region = sys.argv[4]
fname = sys.argv[5]
burn = int(sys.argv[6])
n_blocks_x = n_blocks_y = int(sys.argv[7])
N_nearest = int(sys.argv[8])
grid_lims = getattr(mg, region + '_lims')
nmonths = 288 
start_year = 1985
mask_name = 'landSea-e'
relp=1e-3

hf = tb.openFile(fname)
n_total = len(hf.root.chain0.PyMCsamples)
indices = np.array(np.linspace(burn, n_total, n_jobs+1), dtype=int)
my_start = indices[i]
my_end = indices[i+1]

infile_base = fname.split('/')[-1].replace('.hdf5','')
outfile_name = '/share/scratch/malaria-atlas-project/MAP-outputs/realizations_nb_%i_%s.hdf5'%(n_blocks_x,'_'.join([infile_base, 'iterations', str(i*iter_per_job), str((i+1)*iter_per_job)]))
# print ('create_many_realizations(%i, %i, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, n_blocks_x, n_blocks_y, %s, relp, mask_fname, n_in_trace=%i)'%(my_start, iter_per_job, outfile_name, my_end))

print 'i: %i'%i
print 'iter_per_job: %i'%i
print 'n_jobs: %i'%i
print 'region: %s'%s
print 'fname: %s'%s
print 'burn: %i'%i
print 'nmonths: %i'%n_months
print 'start_year: %i'%start_year
print 'mask_name: %s'%mask_name
print 'relp: %f'%relp
print 'grid_lims: %s'%str(grid_lims)
print 'n_blocks_x: %i'%n_blocks_x
print 'n_blocks_y: %i'%n_blocks_y
print 'N_nearest: %i'%N_nearest

create_many_realizations(my_start, iter_per_job, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, n_blocks_x, n_blocks_y, outfile_name, N_nearest, relp, mask_name, n_in_trace = my_end)
