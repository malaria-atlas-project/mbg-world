# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import numpy as np
import pymc as pm
from pymc.gp.incomplete_chol import ichol_full
from numpy import float32
import tables as tb
import time
import scipy

pm.__PyMCThreadPool__.setNumWorkers(0)

__all__ = ['preprocess', 'krige_month', 'ndmeshgrid']

def ndmeshgrid(grids, hnode=None):
    """
    Converts a list of (start, stop, n) tuples to an 'n-dimensional meshgrid'.
    In two dimensions, this would be:
    
        x = linspace(*grids[0])
        y = linspace(*grids[1])
        x,y = meshgrid(x,y)
        z = concatenate(x,y,axis=-1)
    
    or something like that. Also returns the number of locations in each direction
    as a list.
    """
    ndim = len(grids)
    grids = np.asarray(grids)
    ns = grids[:,2]
    axes = [np.linspace(*grid) for grid in grids]
    if hnode is None:
        x = np.empty(list(ns)+[ndim])
        for index in np.ndindex(*ns):
            x[index+(None,)] = [axes[i][index[i]] for i in xrange(ndim)]
        return np.atleast_2d(x.squeeze()), ns            
    else:
        for index in np.ndindex(*ns):
            hnode[index] = [axes[i][index[i]] for i in xrange(ndim)]
        return ns


def preprocess(C, data_locs, grids, x, n_blocks_x, n_blocks_y, tdata, pdata, relp, mean_ondata): 

    xbi = np.asarray(np.linspace(0,grids[0][2],n_blocks_x+1),dtype=int)
    ybi = np.asarray(np.linspace(0,grids[1][2],n_blocks_y+1),dtype=int)

    dev = (tdata-mean_ondata-pdata)
    
    print 'Max dev:', np.max(np.abs(dev))

    # Figure out which data locations are relevant to the different prediction blocks
    cutoff = C.params['amp']*relp
    scale = C.params['scale']
    inc = C.params['inc']
    ecc = C.params['ecc']
    eff_spat_scale = scale/np.sqrt(2)
    rel_data_ind = np.empty((n_blocks_x, n_blocks_y), dtype=object)
    
    C_eval = C(data_locs,data_locs)
    
    U = np.asarray(np.linalg.cholesky(C_eval).T, order='F')
    dl_posdef = data_locs
    dev_posdef = dev
    # U, n_posdef, pivots = ichol_full(c=C_eval, reltol=relp)
    # U = U[:n_posdef, :n_posdef]
    # 
    # dl_posdef = data_locs[pivots[:n_posdef]]
    # dev_posdef = dev[pivots[:n_posdef]]

    # Backsolve data-data covariance against dev
    pm.gp.trisolve(U, dev_posdef, uplo='U', transa='T', inplace=True)
    pm.gp.trisolve(U, dev_posdef, uplo='U', transa='N', inplace=True)

    return dev_posdef, xbi, ybi, dl_posdef
        

def krige_month(C, i, dl_posdef, grid_shape, n_blocks_x, n_blocks_y, xbi, ybi, x, dev_posdef, row, mask):
        
    x_index_start = 0

    for j in xrange(n_blocks_x):
        x_block_size = (xbi[j+1]-xbi[j])
        
        for k in xrange(n_blocks_y):
            
            this_x = x[xbi[j]:xbi[j+1], ybi[k]:ybi[k+1]]
            this_mask = mask[xbi[j]:xbi[j+1],ybi[k]:ybi[k+1]]
            y_block_size = ybi[k+1]-ybi[k]            
            
            # Check if this block contains any land area, otherwise leave the covariance block as zero.
            if np.sum(this_mask) > 0:   

                this_C_V = C(dl_posdef, this_x)
                this_mask = mask[xbi[j]:xbi[j+1],ybi[k]:ybi[k+1]]             

                row[xbi[j]:xbi[j+1],ybi[k]:ybi[k+1]] = scipy.linalg.blas.fblas.sgemv(1., this_C_V.T, dev_posdef).reshape((x_block_size, y_block_size))

    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()
