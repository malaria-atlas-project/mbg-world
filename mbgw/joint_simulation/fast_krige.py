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


        

