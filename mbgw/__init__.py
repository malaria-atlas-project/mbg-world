# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

#from testmbgw import test
try:
    from testmbgw import test
except ImportError:
    pass

from model import *
import pymc as pm
import numpy as np
import os
from copy import copy
from correction_factors import age_corr_likelihoods, age_corr_factors, two_ten_factors
from scipy import interpolate as interp
from st_cov_fun import *
import time
import auxiliary_data
# import MAPData
import gc
# from get_covariates import extract_environment_to_hdf5
from tables import ObjectAtom
from generic_mbg import *


# Extra stuff for predictive ops.
n_facs = 1000
# postproc = invlogit
def postproc(x=None, lo_age=None, up_age=None, two_ten_facs=two_ten_factors(n_facs)):
    if lo_age is not None:
        facs = age_corr_factors(lo_age, up_age, n_facs)
        def postproc_(x, lo_age=lo_age, up_age=up_age, facs=facs):
            return pm.flib.invlogit(x) * facs[:,np.random.randint(n_facs)]
        return postproc_
    else:
        return invlogit(x) * two_ten_facs[np.random.randint(n_facs)]
    
metadata_keys = ['ti','fi','ui','with_stukel','chunk','disttol','ttol']
f_name = 'eps_p_f'
nugget_name = 'V'
f_has_nugget = True
x_name = 'data_mesh'
diag_safe = True
non_cov_columns = {'lo_age': 'int', 'up_age': 'int', 'pos': 'float', 'neg': 'float'}

bins = np.array([0,.001,.01,.05,.1,.2,.4,1])

def binfn(arr, bins=bins):
    out = np.digitize(arr, bins)
    return out

bin_reduce = histogram_reduce(bins,binfn)

def bin_finalize(products, n, bins=bins, bin_reduce=bin_reduce):
    out = {}
    for i in xrange(len(bins)-1):
        out['p-class-%i-%i'%(bins[i]*100,bins[i+1]*100)] = products[bin_reduce][:,i+1].astype('float')/n
    out['most-likely-class'] = np.argmax(products[bin_reduce], axis=1)
    out['p-most-likely-class'] = np.max(products[bin_reduce], axis=1).astype('float') / n
    return out
        
extra_reduce_fns = [bin_reduce]    
extra_finalize = bin_finalize