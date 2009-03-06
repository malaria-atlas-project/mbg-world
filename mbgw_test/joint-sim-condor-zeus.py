import numpy as np
import pymc as pm
import mbgw
from mbgw.joint_simulation import *
import tables as tb
from mbgw.master_grid import *
import sys

print 'Imports done'

# Establish blocks based on spatial distance only.
n_blocks_x = 20
n_blocks_y = 20

# i = int(sys.argv[1])
# iter_per_job = int(sys.argv[2])
# n_jobs = int(sys.argv[3])

i = 0
iter_per_job = 3
n_jobs = 5

def simulate(fname, grid_lims, nmonths, start_year, burn, mask_name, relp=1e-3):
    """
    n = number of realizations desired
    fname = path to trace file
    grid_lims = [[lon_llc lon_urc] [lat_llc lat_urc]]
    nmonths = number of months to realize
    start_month = starting year
    burn = number of burnin iterations to discard
    outfile_name = name of hdf5 archive
    relp = covariance's relative precision, make it pretty big for speed.
    """
    hf = tb.openFile(fname)
    n_total = len(hf.root.chain0.PyMCsamples)
    indices = np.array(np.linspace(burn, n_total, n_jobs+1), dtype=int)
    my_start = indices[i]
    my_end = indices[i+1]

    infile_base = fname.split('/')[-1].replace('.hdf5','')
    outfile_name = '/share/scratch/zoo-images/MAP-outputs/realizations_%s.hdf5'%('_'.join([infile_base]+[str(j) for j in indices]))
    print ('create_many_realizations(%i, %i, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, n_blocks_x, n_blocks_y, %s, relp, mask_fname, n_in_trace=%i)'%(my_start, iter_per_job, outfile_name, my_end))

    # create_many_realizations(my_start, iter_per_job, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, n_blocks_x, n_blocks_y, outfile_name, relp, mask_name, n_in_trace = my_end)



# Generate simulation
simulate( '../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5', AF_lims, 288, 1985, 2000, 'landSea-e')