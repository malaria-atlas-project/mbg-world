import numpy as np
import pymc as pm
import mbgw
from mbgw.joint_simulation import *
import tables as tb
from mbgw.master_grid import *

# Maximum memory usage in kriging stage
memmax = 1.e8

def simulate(n, fname, grid_lims, nmonths, start_year, burn, outfile_name, mask_name, relp=1e-3):
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
    return create_many_realizations(burn,n, hf.root.chain0, hf.root.metadata, grid_lims, start_year, nmonths, outfile_name, memmax, relp, mask_name)

if __name__ == '__main__':
    try:
        hf.close()
    except:
        pass
    simulate(1, '../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5', AF_lims, 288, 1985, 2000, 'test_sim_%f.hdf5'%memmax, 'landSea-e')

    # import pylab as pl
    # pl.clf()
    # hr = tb.openFile('test_sim.hdf5').root
    # pl.imshow(hr.realizations[0,:,:,0], extent=[hr.axes[0][0], hr.axes[0][1], hr.axes[1][0], hr.axes[1][1]], interpolation='nearest')
    # 
    # pl.colorbar()
